# 0004. 固定依赖

[English](0004-pinned-dependencies.md) | [日本語](0004-pinned-dependencies.ja.md) | **简体中文**

- 状态: Accepted
- 日期: 2026-11-27

## 背景 (Context)
要保证可复现,测试和基准测试在每台机器上、在 CI 中,都必须使用相同版本的库。

## 决定 (Decision)
- C++:依赖通过 CMake 的 `FetchContent` 获取,同时指定发布标签**和** SHA-256 (`cmake/Dependencies.cmake`):Catch2 v3.8.1、Google Benchmark v1.9.4。
- Python:由 uv 管理;依赖用 `uv add` 添加(工具类用 `uv add --dev`),并锁定在已提交的 `uv.lock` 中。不使用 pip。
- 升级固定的版本,需要一份说明理由的新 ADR(或修订用的 ADR)。

## 后果 (Consequences)
- FetchContent 缓存建立之后,构建可以离线复现。
- 版本升级是有意的,并经过评审。
