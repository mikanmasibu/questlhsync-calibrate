#include "calib.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <map>
#include <mutex>

#include "json.h"

static constexpr double kCapture = 12.0;  // s of looking around
static constexpr int kMinSamples = 40;

static Quat Conj(Quat q) { return {q.w, -q.x, -q.y, -q.z}; }

static V3 RotVec(Quat from, Quat to) {  // world rotation taking `from` onto `to`, radians
  Quat d = to * Conj(from);
  if (d.w < 0) d = {-d.w, -d.x, -d.y, -d.z};
  double n = std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
  double ang = 2 * std::atan2(n, std::max(0.0, d.w));
  if (n < 1e-12) return {};
  return V3{d.x, d.y, d.z} * (ang / n);
}

static M3 ExpRot(V3 w) {
  double th = norm(w);
  M3 R;
  if (th < 1e-12) {
    R.m[0][1] = -w.z; R.m[0][2] = w.y;
    R.m[1][0] = w.z;  R.m[1][2] = -w.x;
    R.m[2][0] = -w.y; R.m[2][1] = w.x;
    return R;
  }
  V3 u = w * (1 / th);
  double c = std::cos(th), s = std::sin(th), k = 1 - c;
  R.m[0][0] = c + k * u.x * u.x;     R.m[0][1] = k * u.x * u.y - s * u.z; R.m[0][2] = k * u.x * u.z + s * u.y;
  R.m[1][0] = k * u.y * u.x + s * u.z; R.m[1][1] = c + k * u.y * u.y;     R.m[1][2] = k * u.y * u.z - s * u.x;
  R.m[2][0] = k * u.z * u.x - s * u.y; R.m[2][1] = k * u.z * u.y + s * u.x; R.m[2][2] = c + k * u.z * u.z;
  return R;
}

static bool Solve9(double A[9][9], double b[9], double x[9]) {
  double M[9][10];
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) M[i][j] = A[i][j];
    M[i][9] = b[i];
  }
  for (int c = 0; c < 9; c++) {
    int p = c;
    for (int r = c + 1; r < 9; r++)
      if (std::fabs(M[r][c]) > std::fabs(M[p][c])) p = r;
    if (std::fabs(M[p][c]) < 1e-14) return false;
    if (p != c)
      for (int j = c; j < 10; j++) std::swap(M[p][j], M[c][j]);
    for (int r = 0; r < 9; r++) {
      if (r == c) continue;
      double f = M[r][c] / M[c][c];
      for (int j = c; j < 10; j++) M[r][j] -= f * M[c][j];
    }
  }
  for (int i = 0; i < 9; i++) x[i] = M[i][9] / M[i][i];
  return true;
}

// Smallest and largest eigenvalues of a 3x3 symmetric matrix, by a few Jacobi sweeps.
static void EigRange(const double S[3][3], double &lo, double &hi) {
  double a[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) a[i][j] = S[i][j];
  for (int sweep = 0; sweep < 8; sweep++)
    for (int p = 0; p < 3; p++)
      for (int q = p + 1; q < 3; q++) {
        if (std::fabs(a[p][q]) < 1e-18) continue;
        double th = (a[q][q] - a[p][p]) / (2 * a[p][q]);
        double t = (th >= 0 ? 1 : -1) / (std::fabs(th) + std::sqrt(th * th + 1));
        double c = 1 / std::sqrt(t * t + 1), s = t * c;
        for (int k = 0; k < 3; k++) {
          double akp = a[k][p], akq = a[k][q];
          a[k][p] = c * akp - s * akq;
          a[k][q] = s * akp + c * akq;
        }
        for (int k = 0; k < 3; k++) {
          double apk = a[p][k], aqk = a[q][k];
          a[p][k] = c * apk - s * aqk;
          a[q][k] = s * apk + c * aqk;
        }
      }
  lo = hi = a[0][0];
  for (int i = 1; i < 3; i++) { lo = std::min(lo, a[i][i]); hi = std::max(hi, a[i][i]); }
}

static double Cost(const std::vector<CalSample> &s, const std::vector<std::pair<int, int>> &pairs, const M3 &R, V3 t,
                   V3 c) {
  double ss = 0;
  for (auto &e : s) {
    V3 r = R * e.pd + t - e.ph - ToM3(e.qh) * c;
    ss += dot(r, r);
  }
  for (auto &ij : pairs) {
    V3 wh = RotVec(s[ij.first].qh, s[ij.second].qh);
    V3 wd = RotVec(s[ij.first].qd, s[ij.second].qd);
    V3 r = R * wd - wh;
    ss += 0.25 * dot(r, r);  // a radian of axis error weighs like 0.5 m
  }
  return ss;
}

CalFit SolveCalibration(const std::vector<CalSample> &in) {
  CalFit out;
  if ((int)in.size() < kMinSamples) { out.err = "not enough motion"; return out; }
  std::vector<CalSample> s = in;
  if (s.size() > 240) {
    size_t step = s.size() / 240;
    std::vector<CalSample> t;
    for (size_t i = 0; i < s.size(); i += step) t.push_back(s[i]);
    s.swap(t);
  }
  V3 pmid{};
  for (auto &e : s) pmid += e.ph;
  pmid = pmid * (1.0 / s.size());
  double reach = 0;  // how far the device sat from the headset: a hip tracker is not a head mount
  for (auto &e : s) reach += norm(e.pd - e.ph);
  reach /= s.size();
  if (reach > 0.40) { out.err = "no tracker or controller stayed on the headset"; return out; }

  std::vector<std::pair<int, int>> pairs;
  double info[3][3] = {};
  for (size_t i = 0; i + 1 < s.size(); i++) {
    for (size_t j = i + 1; j < s.size() && j < i + 30; j++) {
      double ang = QuatDeg(s[i].qh, s[j].qh);
      if (ang < 8 || ang > 50) continue;
      V3 wd = RotVec(s[i].qd, s[j].qd);
      // turning the head about wd observes the correction's rotation perpendicular to wd
      double d = dot(wd, wd);
      for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) info[a][b] += (a == b ? d : 0) - wd[a] * wd[b];
      pairs.push_back({(int)i, (int)j});
      break;
    }
  }
  double lo, hi;
  EigRange(info, lo, hi);
  if (pairs.size() < 12 || hi < 0.15) { out.err = "look further: left, right, up and down"; return out; }

  M3 R;
  V3 t{}, c{};
  for (auto &e : s) c += T(ToM3(e.qh)) * (e.pd - e.ph);
  c = c * (1.0 / s.size());
  double lam = 1e-3, prev = Cost(s, pairs, R, t, c);
  for (int it = 0; it < 12; it++) {
    double A[9][9] = {}, b[9] = {};
    auto add = [&](const double J[9], const V3 &r, double w) {
      for (int k = 0; k < 3; k++) {
        for (int i = 0; i < 9; i++) {
          b[i] += w * J[k * 9 + i] * r[k];
          for (int j = 0; j < 9; j++) A[i][j] += w * J[k * 9 + i] * J[k * 9 + j];
        }
      }
    };
    for (auto &e : s) {
      M3 Rh = ToM3(e.qh);
      V3 rp = R * e.pd;
      V3 r = rp + t - e.ph - Rh * c;
      // left perturbation R <- Exp(w) R: d(Rp) = w cross Rp
      double J[27] = {};
      J[0 * 9 + 1] = rp.z;  J[0 * 9 + 2] = -rp.y;
      J[1 * 9 + 0] = -rp.z; J[1 * 9 + 2] = rp.x;
      J[2 * 9 + 0] = rp.y;  J[2 * 9 + 1] = -rp.x;
      J[0 * 9 + 3] = J[1 * 9 + 4] = J[2 * 9 + 5] = 1;
      for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++) J[row * 9 + 6 + col] = -Rh.m[row][col];
      add(J, r, 1);
    }
    for (auto &ij : pairs) {
      V3 wh = RotVec(s[ij.first].qh, s[ij.second].qh);
      V3 wd = RotVec(s[ij.first].qd, s[ij.second].qd);
      V3 rw = R * wd;
      V3 r = rw - wh;
      double J[27] = {};
      J[0 * 9 + 1] = rw.z;  J[0 * 9 + 2] = -rw.y;
      J[1 * 9 + 0] = -rw.z; J[1 * 9 + 2] = rw.x;
      J[2 * 9 + 0] = rw.y;  J[2 * 9 + 1] = -rw.x;
      add(J, r, 0.25);
    }
    for (int i = 0; i < 9; i++) A[i][i] += lam;
    double rhs[9], dx[9];
    for (int i = 0; i < 9; i++) rhs[i] = -b[i];
    if (!Solve9(A, rhs, dx)) { lam *= 10; continue; }
    double step = 0;
    for (int i = 0; i < 9; i++) step += std::fabs(dx[i]);
    M3 Rn = ExpRot(V3{dx[0], dx[1], dx[2]}) * R;
    V3 tn = t + V3{dx[3], dx[4], dx[5]}, cn = c + V3{dx[6], dx[7], dx[8]};
    double now = Cost(s, pairs, Rn, tn, cn);
    if (now <= prev) {
      R = Rn; t = tn; c = cn; prev = now;
      lam = std::max(lam / 3, 1e-8);
      if (step < 1e-7) break;
    } else {
      lam *= 4;
      if (lam > 1e6) break;
    }
  }
  double ss = 0;
  for (auto &e : s) {
    V3 r = R * e.pd + t - e.ph - ToM3(e.qh) * c;
    ss += dot(r, r);
  }
  out.rms = std::sqrt(ss / s.size());
  if (out.rms > 0.015) { out.err = "the tracker moved against the headset; hold it firmly and try again"; return out; }
  out.deg = RotDeg(R);
  V3 shift = R * pmid + t - pmid;
  out.cm = norm(shift) * 100;
  if (out.deg > 8 || out.cm > 8) { out.err = "the correction is too large to trust; was the tracker fixed to the headset?"; return out; }
  out.q = ToQuat(R);
  out.t = t;
  out.ok = true;
  return out;
}

// ---------------------------------------------------------------- capture
Calibrator::Calibrator(std::string path, LogFn log) : path_(std::move(path)), log_(std::move(log)) {}

void Calibrator::Write(bool on, Quat q, V3 t, double deg, double cm) {
  seq_.fetch_add(1, std::memory_order_acq_rel);
  std::atomic_thread_fence(std::memory_order_release);
  on_ = on;
  q_ = q;
  t_ = t;
  deg_ = deg;
  cm_ = cm;
  phase_ = on ? 2 : 0;
  std::atomic_thread_fence(std::memory_order_release);
  seq_.fetch_add(1, std::memory_order_acq_rel);
}

bool Calibrator::Read(Quat &q, V3 &t) const {
  uint32_t s = 0;
  Quat cq;
  V3 ct;
  int st = 0;
  for (int i = 0; i < 8; i++) {
    s = seq_.load(std::memory_order_acquire);
    if (s & 1) continue;
    cq = q_;
    ct = t_;
    st = on_;
    std::atomic_thread_fence(std::memory_order_acquire);
    if (seq_.load(std::memory_order_acquire) == s) break;
  }
  if (!st) return false;
  q = cq;
  t = ct;
  return true;
}

void Calibrator::Save() const {
  FILE *f = fopen(path_.c_str(), "w");
  if (!f) return;
  fprintf(f, "{\n \"q\": [%.9f, %.9f, %.9f, %.9f],\n \"t\": [%.6f, %.6f, %.6f],\n \"deg\": %.4f,\n \"cm\": %.3f\n}\n",
          q_.w, q_.x, q_.y, q_.z, t_.x, t_.y, t_.z, deg_, cm_);
  fclose(f);
}

void Calibrator::Load() {
  FILE *f = fopen(path_.c_str(), "rb");
  if (!f) return;
  std::string text;
  char buf[512];
  size_t n;
  while ((n = fread(buf, 1, sizeof buf, f)) > 0) text.append(buf, n);
  fclose(f);
  JVal d;
  if (!JParse(text, d)) return;
  const JVal *q = d.get("q"), *tv = d.get("t");
  if (!q || !tv) return;
  auto qn = q->nums(), tn = tv->nums();
  if (qn.size() < 4 || tn.size() < 3) return;
  Quat cq{qn[0], qn[1], qn[2], qn[3]};
  double nn = std::sqrt(cq.w * cq.w + cq.x * cq.x + cq.y * cq.y + cq.z * cq.z);
  if (!(nn > 0.5 && nn < 1.5)) return;
  cq = {cq.w / nn, cq.x / nn, cq.y / nn, cq.z / nn};
  V3 ct{tn[0], tn[1], tn[2]};
  double deg = d.get("deg") ? d.get("deg")->num() : RotDeg(ToM3(cq));
  double cm = d.get("cm") ? d.get("cm")->num() : norm(ct) * 100;
  Write(true, cq, ct, deg, cm);
  char msg[160];
  snprintf(msg, sizeof msg, "calibration from last session: %.2f deg, %.1f cm at the head", deg, cm);
  log_(msg);
}

bool Calibrator::Start(double now) {
  if (sampling_.load(std::memory_order_acquire)) return false;
  {
    std::lock_guard<std::mutex> g(mu_);
    dev_.clear();
    hmd_ = {};
  }
  t0_ = now;
  pct_ = 0;
  phase_ = 1;
  sampling_.store(true, std::memory_order_release);
  log_("calibration: hold a tracker or controller firmly on the headset and look left, right, up and down");
  return true;
}

void Calibrator::Cancel() {
  if (!sampling_.exchange(false)) return;
  {
    std::lock_guard<std::mutex> g(mu_);
    dev_.clear();
  }
  pct_ = 0;
  phase_ = on_ ? 2 : 0;
  log_("calibration cancelled");
}

void Calibrator::Clear() {
  if (sampling_.exchange(false)) {
    std::lock_guard<std::mutex> g(mu_);
    dev_.clear();
  }
  Write(false, Quat{1, 0, 0, 0}, {}, 0, 0);
  remove(path_.c_str());
  log_("calibration cleared");
}

bool Calibrator::Tick(double now) {
  if (!sampling_.load(std::memory_order_acquire)) return false;
  pct_ = (int)std::min(100.0, (now - t0_) / kCapture * 100);
  if (now - t0_ < kCapture) return false;
  sampling_.store(false, std::memory_order_release);
  Finish();
  return true;
}

void Calibrator::Finish() {
  std::map<int, Dev> dev;
  {
    std::lock_guard<std::mutex> g(mu_);
    dev.swap(dev_);
  }
  const Dev *best = nullptr;
  double best_d = 1e9;
  for (auto &kv : dev) {
    if ((int)kv.second.s.size() < kMinSamples) continue;
    double d = 0;
    for (auto &e : kv.second.s) d += norm(e.pd - e.ph);
    d /= kv.second.s.size();
    if (d < best_d) { best_d = d; best = &kv.second; }
  }
  pct_ = 0;
  auto restore = [&] {
    Quat q; V3 t;
    phase_ = Read(q, t) ? 2 : 0;
  };
  if (!best) {
    restore();
    log_("calibration: no tracker or controller stayed on the headset");
    return;
  }
  CalFit fit = SolveCalibration(best->s);
  if (!fit.ok) {
    restore();
    log_(std::string("calibration: ") + fit.err);
    return;
  }
  if (fit.deg < 0.05 && fit.cm < 0.15) {
    restore();
    log_("calibration: already within 0.05 deg and 1.5 mm, left as it is");
    return;
  }
  Quat oq{1, 0, 0, 0};
  V3 ot{};
  bool had = Read(oq, ot);
  M3 R = ToM3(fit.q) * (had ? ToM3(oq) : M3());
  V3 t = ToM3(fit.q) * ot + fit.t;
  V3 pmid{};
  for (auto &e : best->s) pmid += e.ph;
  pmid = pmid * (1.0 / best->s.size());
  double deg = RotDeg(R), cm = norm(R * pmid + t - pmid) * 100;
  Write(true, ToQuat(R), t, deg, cm);
  Save();
  char msg[200];
  snprintf(msg, sizeof msg, "calibration: %.2f deg, %.1f cm at the head (the mount held to %.0f mm)", fit.deg, fit.cm,
           fit.rms * 1000);
  log_(msg);
}

void Calibrator::OnHmd(double t, Quat q, V3 p) {
  if (!sampling_.load(std::memory_order_acquire)) return;
  std::lock_guard<std::mutex> g(mu_);
  hmd_ = {t, q, p, true};
}

void Calibrator::OnDevice(int id, bool reference, double t, Quat q, V3 p) {
  if (reference || !sampling_.load(std::memory_order_acquire)) return;
  if (!std::isfinite(p.x + p.y + p.z)) return;
  std::lock_guard<std::mutex> g(mu_);
  if (!hmd_.has || std::fabs(t - hmd_.t) > 0.04) return;
  Dev &d = dev_[id];
  if (d.s.size() >= 600) return;
  if (!d.s.empty() && t - d.t < 0.03) return;
  if (!d.s.empty() && QuatDeg(d.q, hmd_.q) < 2.0 && norm(p - d.p) < 0.01) return;
  d.t = t;
  d.q = hmd_.q;
  d.p = p;
  d.s.push_back({hmd_.q, q, hmd_.p, p});
}

#ifdef CALIB_SELFTEST
#include <cstdio>
#include <cstdlib>
static int Fail(const char *m) { fprintf(stderr, "%s\n", m); return 1; }
int main() {
  V3 ctrue{0.02, -0.08, 0.04};
  M3 Rm = ExpRot(V3{0.2, -0.4, 0.1});
  M3 Rerr = ExpRot(V3{1.5 / kDeg, -0.4 / kDeg, 0.7 / kDeg});  // tilt plus a little yaw
  V3 terr{0.012, -0.008, 0.015};
  std::vector<CalSample> s;
  for (int iy = -2; iy <= 2; iy++)
    for (int ip = -1; ip <= 1; ip++)
      for (int ir = -1; ir <= 1; ir++) {
        M3 Rh = ExpRot(V3{ip * 18.0 / kDeg, iy * 25.0 / kDeg, ir * 12.0 / kDeg});
        V3 ph{0.4 + iy * 0.02, 1.6 + ip * 0.01, 1.1 + ir * 0.02};
        V3 ptr = ph + Rh * ctrue;
        M3 Rd = Rh * Rm;
        CalSample e;
        e.qh = ToQuat(Rh);
        e.ph = ph;
        e.pd = T(Rerr) * (ptr - terr);
        e.qd = ToQuat(T(Rerr) * Rd);
        s.push_back(e);
      }
  CalFit f = SolveCalibration(s);
  if (!f.ok) return Fail(f.err);
  double dang = RotDeg(T(ToM3(f.q)) * Rerr);
  double dtr = norm(ToM3(f.q) * terr + f.t - terr);  // applied to a point: checked below
  // Rfit * p_obs + tfit should equal p_true, and p_obs = Rerr^T * (p_true - terr),
  // so Rfit * Rerr^T * (p - terr) + tfit = p for every p. Thus Rfit = Rerr and tfit = terr.
  double rerr = RotDeg(T(Rerr) * ToM3(f.q));
  double terr_err = norm(f.t - terr);
  printf("recovered %.3f deg (err %.4f deg) t err %.3f mm  rms %.2f mm  n %d\n", f.deg, rerr, terr_err * 1000,
         f.rms * 1000, (int)s.size());
  if (rerr > 0.05 || terr_err > 0.001) return Fail("solver missed the correction");
  (void)dang; (void)dtr;
  // identity: no error, must report already-small and ok
  for (auto &e : s) {
    M3 Rh = ToM3(e.qh);
    e.pd = e.ph + Rh * ctrue;
    e.qd = ToQuat(Rh * Rm);
  }
  CalFit id = SolveCalibration(s);
  if (!id.ok) return Fail(id.err);
  printf("identity %.4f deg %.3f mm\n", id.deg, id.cm * 10);
  if (id.deg > 0.02 || id.cm > 0.05) return Fail("identity was not identity");
  puts("ok");
  return 0;
}
#endif
