// Minimal strict JSON (RFC 8259) for the benchmark runner: read --env-json, write records.

#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>

#include "sssp/bench_runner.hpp"

namespace sssp::bench {

namespace {

constexpr int kMaxDepth = 256;

bool is_digit(char c) noexcept {
  return c >= '0' && c <= '9';
}

/// End of the JSON number starting at `pos`, or npos if there is no valid number there.
std::size_t scan_number(std::string_view s, std::size_t pos) noexcept {
  auto at = [&](std::size_t i) { return i < s.size() ? s[i] : '\0'; };
  if (at(pos) == '-') ++pos;
  if (at(pos) == '0') {
    ++pos;
  } else if (is_digit(at(pos))) {
    while (is_digit(at(pos))) ++pos;
  } else {
    return std::string_view::npos;
  }
  if (at(pos) == '.') {
    ++pos;
    if (!is_digit(at(pos))) return std::string_view::npos;
    while (is_digit(at(pos))) ++pos;
  }
  if (at(pos) == 'e' || at(pos) == 'E') {
    ++pos;
    if (at(pos) == '+' || at(pos) == '-') ++pos;
    if (!is_digit(at(pos))) return std::string_view::npos;
    while (is_digit(at(pos))) ++pos;
  }
  if (is_digit(at(pos))) return std::string_view::npos;  // leading zero, e.g. 01
  return pos;
}

class Parser {
 public:
  explicit Parser(std::string_view s) : s_(s) {}

  Json document() {
    Json v = value(0);
    skip_ws();
    if (pos_ != s_.size()) fail("trailing characters after the JSON value");
    return v;
  }

 private:
  [[noreturn]] void fail(const std::string& what) const {
    throw JsonError("invalid JSON at offset " + std::to_string(pos_) + ": " + what);
  }

  void skip_ws() noexcept {
    while (pos_ < s_.size() &&
           (s_[pos_] == ' ' || s_[pos_] == '\t' || s_[pos_] == '\n' || s_[pos_] == '\r')) {
      ++pos_;
    }
  }

  [[nodiscard]] char peek() const noexcept { return pos_ < s_.size() ? s_[pos_] : '\0'; }

  void expect(char c) {
    if (peek() != c) fail(std::string("expected '") + c + "'");
    ++pos_;
  }

  void literal(std::string_view word) {
    if (s_.substr(pos_, word.size()) != word) fail("invalid literal");
    pos_ += word.size();
  }

  Json value(int depth) {
    if (depth > kMaxDepth) fail("nesting too deep");
    skip_ws();
    const char c = peek();
    switch (c) {
      case '{':
        return object(depth);
      case '[':
        return array(depth);
      case '"':
        return Json::string(string());
      case 't':
        literal("true");
        return Json::boolean(true);
      case 'f':
        literal("false");
        return Json::boolean(false);
      case 'n':
        literal("null");
        return Json::null();
      default:
        if (c == '-' || is_digit(c)) return number();
        fail("unexpected character");
    }
  }

  Json object(int depth) {
    expect('{');
    Json obj = Json::object();
    skip_ws();
    if (peek() == '}') {
      ++pos_;
      return obj;
    }
    for (;;) {
      skip_ws();
      if (peek() != '"') fail("expected a string key");
      std::string key = string();
      skip_ws();
      expect(':');
      Json v = value(depth + 1);
      obj.set(key, std::move(v));  // duplicate keys: last value wins (Python json semantics)
      skip_ws();
      if (peek() == ',') {
        ++pos_;
        continue;
      }
      expect('}');
      return obj;
    }
  }

  Json array(int depth) {
    expect('[');
    Json arr = Json::array();
    skip_ws();
    if (peek() == ']') {
      ++pos_;
      return arr;
    }
    for (;;) {
      arr.push_back(value(depth + 1));
      skip_ws();
      if (peek() == ',') {
        ++pos_;
        continue;
      }
      expect(']');
      return arr;
    }
  }

  Json number() {
    const std::size_t end = scan_number(s_, pos_);
    if (end == std::string_view::npos) fail("invalid number");
    const std::string_view token = s_.substr(pos_, end - pos_);
    pos_ = end;
    return Json::raw_number(token);
  }

  unsigned hex4() {
    if (pos_ + 4 > s_.size()) fail("truncated \\u escape");
    unsigned v = 0;
    for (int i = 0; i < 4; ++i) {
      const char c = s_[pos_++];
      v <<= 4U;
      if (c >= '0' && c <= '9') {
        v |= static_cast<unsigned>(c - '0');
      } else if (c >= 'a' && c <= 'f') {
        v |= static_cast<unsigned>(c - 'a' + 10);
      } else if (c >= 'A' && c <= 'F') {
        v |= static_cast<unsigned>(c - 'A' + 10);
      } else {
        fail("invalid \\u escape");
      }
    }
    return v;
  }

  static void append_utf8(std::string& out, unsigned cp) {
    if (cp < 0x80) {
      out += static_cast<char>(cp);
    } else if (cp < 0x800) {
      out += static_cast<char>(0xC0U | (cp >> 6U));
      out += static_cast<char>(0x80U | (cp & 0x3FU));
    } else if (cp < 0x10000) {
      out += static_cast<char>(0xE0U | (cp >> 12U));
      out += static_cast<char>(0x80U | ((cp >> 6U) & 0x3FU));
      out += static_cast<char>(0x80U | (cp & 0x3FU));
    } else {
      out += static_cast<char>(0xF0U | (cp >> 18U));
      out += static_cast<char>(0x80U | ((cp >> 12U) & 0x3FU));
      out += static_cast<char>(0x80U | ((cp >> 6U) & 0x3FU));
      out += static_cast<char>(0x80U | (cp & 0x3FU));
    }
  }

  std::string string() {
    expect('"');
    std::string out;
    for (;;) {
      if (pos_ >= s_.size()) fail("unterminated string");
      const char c = s_[pos_++];
      if (c == '"') return out;
      if (static_cast<unsigned char>(c) < 0x20) fail("control character in string");
      if (c != '\\') {
        out += c;
        continue;
      }
      if (pos_ >= s_.size()) fail("unterminated escape");
      const char e = s_[pos_++];
      switch (e) {
        case '"':
          out += '"';
          break;
        case '\\':
          out += '\\';
          break;
        case '/':
          out += '/';
          break;
        case 'b':
          out += '\b';
          break;
        case 'f':
          out += '\f';
          break;
        case 'n':
          out += '\n';
          break;
        case 'r':
          out += '\r';
          break;
        case 't':
          out += '\t';
          break;
        case 'u': {
          unsigned cp = hex4();
          if (cp >= 0xD800 && cp <= 0xDBFF) {
            if (s_.substr(pos_, 2) != "\\u") fail("unpaired surrogate");
            pos_ += 2;
            const unsigned lo = hex4();
            if (lo < 0xDC00 || lo > 0xDFFF) fail("unpaired surrogate");
            cp = 0x10000 + ((cp - 0xD800) << 10U) + (lo - 0xDC00);
          } else if (cp >= 0xDC00 && cp <= 0xDFFF) {
            fail("unpaired surrogate");
          }
          append_utf8(out, cp);
          break;
        }
        default:
          fail("invalid escape");
      }
    }
  }

  std::string_view s_;
  std::size_t pos_ = 0;
};

template <class T>
std::string to_text(T v) {
  std::array<char, 64> buf{};
  const auto [ptr, ec] = std::to_chars(buf.data(), buf.data() + buf.size(), v);
  if (ec != std::errc{}) throw JsonError("number formatting failed");
  return {buf.data(), ptr};
}

void dump_string(std::string& out, const std::string& s) {
  static constexpr std::string_view kHex = "0123456789abcdef";
  out += '"';
  for (const char c : s) {
    const auto u = static_cast<unsigned char>(c);
    switch (c) {
      case '"':
        out += "\\\"";
        break;
      case '\\':
        out += "\\\\";
        break;
      case '\b':
        out += "\\b";
        break;
      case '\f':
        out += "\\f";
        break;
      case '\n':
        out += "\\n";
        break;
      case '\r':
        out += "\\r";
        break;
      case '\t':
        out += "\\t";
        break;
      default:
        if (u < 0x20) {
          out += "\\u00";
          out += kHex[u >> 4U];
          out += kHex[u & 0xFU];
        } else {
          out += c;
        }
    }
  }
  out += '"';
}

}  // namespace

Json Json::boolean(bool b) {
  Json j;
  j.kind_ = Kind::kBool;
  j.bool_ = b;
  return j;
}

Json Json::integer(std::uint64_t v) {
  return raw_number(to_text(v));
}

Json Json::integer(std::int64_t v) {
  return raw_number(to_text(v));
}

Json Json::real(double v) {
  if (!std::isfinite(v)) throw JsonError("JSON numbers must be finite");
  return raw_number(to_text(v));
}

Json Json::raw_number(std::string_view token) {
  if (token.empty() || scan_number(token, 0) != token.size()) {
    throw JsonError("not a JSON number: '" + std::string(token) + "'");
  }
  Json j;
  j.kind_ = Kind::kNumber;
  j.text_ = std::string(token);
  return j;
}

Json Json::string(std::string s) {
  Json j;
  j.kind_ = Kind::kString;
  j.text_ = std::move(s);
  return j;
}

Json Json::array() {
  Json j;
  j.kind_ = Kind::kArray;
  return j;
}

Json Json::object() {
  Json j;
  j.kind_ = Kind::kObject;
  return j;
}

Json Json::parse(std::string_view text) {
  return Parser(text).document();
}

bool Json::is_integer() const noexcept {
  return kind_ == Kind::kNumber && text_.find_first_of(".eE") == std::string::npos;
}

bool Json::as_bool() const {
  if (kind_ != Kind::kBool) throw JsonError("not a boolean");
  return bool_;
}

const std::string& Json::as_string() const {
  if (kind_ != Kind::kString) throw JsonError("not a string");
  return text_;
}

const std::string& Json::number_text() const {
  if (kind_ != Kind::kNumber) throw JsonError("not a number");
  return text_;
}

const Json::Array& Json::elements() const {
  if (kind_ != Kind::kArray) throw JsonError("not an array");
  return array_;
}

const Json::Object& Json::items() const {
  if (kind_ != Kind::kObject) throw JsonError("not an object");
  return object_;
}

const Json* Json::find(std::string_view key) const {
  for (const auto& [k, v] : items()) {
    if (k == key) return &v;
  }
  return nullptr;
}

Json* Json::find_mut(std::string_view key) {
  if (kind_ != Kind::kObject) throw JsonError("not an object");
  for (auto& [k, v] : object_) {
    if (k == key) return &v;
  }
  return nullptr;
}

void Json::set(std::string_view key, Json value) {
  if (Json* existing = find_mut(key)) {
    *existing = std::move(value);
    return;
  }
  object_.emplace_back(std::string(key), std::move(value));
}

void Json::push_back(Json value) {
  if (kind_ != Kind::kArray) throw JsonError("not an array");
  array_.push_back(std::move(value));
}

std::string Json::dump() const {
  std::string out;
  dump_to(out);
  return out;
}

void Json::dump_to(std::string& out) const {
  switch (kind_) {
    case Kind::kNull:
      out += "null";
      break;
    case Kind::kBool:
      out += bool_ ? "true" : "false";
      break;
    case Kind::kNumber:
      out += text_;
      break;
    case Kind::kString:
      dump_string(out, text_);
      break;
    case Kind::kArray: {
      out += '[';
      bool first = true;
      for (const auto& v : array_) {
        if (!first) out += ',';
        first = false;
        v.dump_to(out);
      }
      out += ']';
      break;
    }
    case Kind::kObject: {
      out += '{';
      bool first = true;
      for (const auto& [k, v] : object_) {
        if (!first) out += ',';
        first = false;
        dump_string(out, k);
        out += ':';
        v.dump_to(out);
      }
      out += '}';
      break;
    }
  }
}

}  // namespace sssp::bench
