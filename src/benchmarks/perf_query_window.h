#pragma once

#include <fcntl.h>
#include <unistd.h>

#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string_view>

namespace pixie::benchmarks {

// Optional query-only perf window for native benchmark executables. perf starts
// disabled and acknowledges each command through the supplied FIFO. Warmup must
// be disabled by the registering harness; cumulative queries/phases let runners
// reject unexpected calibration, warmup, or extra selected rows.
class PerfQueryWindow {
 public:
  static bool requested() {
    return std::getenv("PIXIE_RANK_SELECT_PERF_CTL") != nullptr;
  }

  static PerfQueryWindow& instance() {
    static PerfQueryWindow window;
    return window;
  }

  void start() {
    if (ctl_ >= 0) {
      ++phases_;
      command("enable\n");
    }
  }

  void stop(std::uint64_t queries) {
    if (ctl_ >= 0) {
      command("disable\n");
      queries_ += queries;
    }
  }

  bool active() const { return ctl_ >= 0; }
  std::uint64_t queries() const { return queries_; }
  std::uint64_t phases() const { return phases_; }

 private:
  PerfQueryWindow() {
    const auto* ctl = std::getenv("PIXIE_RANK_SELECT_PERF_CTL");
    if (ctl == nullptr) {
      return;
    }
    const auto* ack = std::getenv("PIXIE_RANK_SELECT_PERF_ACK");
    if (ack == nullptr) {
      throw std::runtime_error("Missing PIXIE_RANK_SELECT_PERF_ACK");
    }
    ctl_ = ::open(ctl, O_WRONLY);
    if (ctl_ < 0) {
      throw std::runtime_error("Cannot open perf control FIFO");
    }
    ack_ = ::open(ack, O_RDONLY);
    if (ack_ < 0) {
      ::close(ctl_);
      throw std::runtime_error("Cannot open perf acknowledgement FIFO");
    }
  }

  ~PerfQueryWindow() {
    if (ctl_ >= 0) {
      ::close(ctl_);
      ::close(ack_);
    }
  }

  void command(std::string_view text) {
    std::size_t done = 0;
    while (done < text.size()) {
      const auto bytes = ::write(ctl_, text.data() + done, text.size() - done);
      if (bytes < 0 && errno == EINTR) {
        continue;
      }
      if (bytes <= 0) {
        throw std::runtime_error("Cannot send perf command");
      }
      done += bytes;
    }
    constexpr char expected[] = "ack\n";
    char response[sizeof(expected)];
    done = 0;
    while (done < sizeof(response)) {
      const auto bytes = ::read(ack_, response + done, sizeof(response) - done);
      if (bytes < 0 && errno == EINTR) {
        continue;
      }
      if (bytes <= 0) {
        throw std::runtime_error("Missing perf acknowledgement");
      }
      done += bytes;
    }
    // perf writes sizeof("ack\n"), including the trailing NUL: five bytes.
    if (std::string_view(response, sizeof(response)) !=
        std::string_view(expected, sizeof(expected))) {
      throw std::runtime_error("Invalid perf acknowledgement");
    }
  }

  int ctl_ = -1;
  int ack_ = -1;
  std::uint64_t queries_ = 0;
  std::uint64_t phases_ = 0;
};

}  // namespace pixie::benchmarks
