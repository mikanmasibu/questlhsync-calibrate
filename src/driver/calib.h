// One-time calibration of the residual between the headset and lighthouse space.
// The camera fit is yaw and translation. A tracker or controller held on the headset, while the head looks
// left, right, up and down, also shows the tilt and the small offset the cameras leave. The mount itself is
// solved for and not applied: only the rigid correction of the lighthouse space is, on top of the live fit.
#pragma once
#include <atomic>
#include <functional>
#include <map>
#include <mutex>
#include <string>
#include <vector>

#include "mathx.h"

struct CalSample { Quat qh, qd; V3 ph, pd; };  // headset and lighthouse device, both in the aligned world

struct CalFit {
  bool ok = false;
  Quat q{1, 0, 0, 0};
  V3 t;
  double deg = 0, cm = 0, rms = 0;  // rotation, shift of the head's position, mount residual (m)
  const char *err = "";
};

CalFit SolveCalibration(const std::vector<CalSample> &samples);

class Calibrator {
 public:
  using LogFn = std::function<void(const std::string &)>;
  Calibrator(std::string path, LogFn log);
  void Load();
  bool sampling() const { return sampling_.load(std::memory_order_acquire); }
  bool Start(double now);  // begins a capture; caller freezes the live fit first
  void Cancel();
  void Clear();
  bool Tick(double now);  // true when a capture just finished (caller resumes the live fit)
  void OnHmd(double t, Quat q, V3 p);
  void OnDevice(int id, bool reference, double t, Quat q, V3 p);
  bool Read(Quat &q, V3 &t) const;  // false: no correction
  int state() const { return phase_; }  // 0 none, 1 sampling, 2 applied
  int pct() const { return pct_; }
  double deg() const { return deg_; }
  double cm() const { return cm_; }

 private:
  std::string path_;
  LogFn log_;
  std::atomic<bool> sampling_{false};
  std::atomic<uint32_t> seq_{0};
  int on_ = 0;  // a correction is stored; the hook reads it through the seqlock
  Quat q_{1, 0, 0, 0};
  V3 t_{};
  int phase_ = 0, pct_ = 0;
  double deg_ = 0, cm_ = 0, t0_ = 0;
  struct Last { double t; Quat q; V3 p; bool has = false; };
  std::mutex mu_;
  Last hmd_;
  struct Dev {
    std::vector<CalSample> s;
    double t = -1e18;
    Quat q{};
    V3 p{};
  };
  std::map<int, Dev> dev_;
  void Write(bool on, Quat q, V3 t, double deg, double cm);
  void Save() const;
  void Finish();
};
