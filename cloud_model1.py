# === Daisyworld + normalized clouds: up/down hysteresis sweep vs Luminosity ===
import numpy as np, math, matplotlib.pyplot as plt

# --- Constants (keep consistent with your earlier code) ---
SIGMA     = 5.670367e-8
S0_TOA    = 1361.0
S0_MODEL  = 917.0
SCALE     = S0_MODEL / (S0_TOA/4.0)

ALB_W, ALB_B, ALB_G = 0.75, 0.25, 0.50
T_MIN, T_MAX = 278.15, 313.15
R      = 0.30     # a bit stronger coupling (was 0.25)
GAMMA  = 0.30     # a bit lower mortality (was 0.25)
DT     = 1.0

# --- Geometry helpers (daily-mean TOA) ---
def deg2rad(d): return d*math.pi/180.0
def declination(day): return deg2rad(23.44)*math.sin(2*math.pi*(day-81)/365.0)
def eccentricity_correction(day): return 1.0 + 0.033*math.cos(2*math.pi*day/365.0)
def daily_mean_toa(lat_deg, day):
    phi = deg2rad(lat_deg); delta = declination(day); E0 = eccentricity_correction(day)
    arg = max(-1.0, min(1.0, -math.tan(phi)*math.tan(delta)))
    omega_s = math.acos(arg)
    return (S0_TOA/math.pi)*E0*(omega_s*math.sin(phi)*math.sin(delta)
                                + math.cos(phi)*math.cos(delta)*math.sin(omega_s))

# Planetary mean (area-weighted cosφ)
lat_edges = np.linspace(-90,90,181)
lat_centers = 0.5*(lat_edges[:-1]+lat_edges[1:])
weights = np.cos(np.deg2rad(lat_centers)); weights /= weights.sum()
def planetary_daily_mean(day):
    vals = np.array([daily_mean_toa(lat, day) for lat in lat_centers])
    return float(np.sum(vals*weights))   # ~ S0_TOA*E0/4

# --- Clouds: monthly K_T normalized to mean 1 ---
# --- Clouds: monthly cycle with global mean cloud fraction ~0.66 ---
# Global mean cloud fraction has a small seasonal cycle (~±0.02–0.03 around 0.66).
CLOUD_FRAC_MONTH = np.array([
    0.66,  # Jan
    0.65,  # Feb
    0.67,  # Mar
    0.66,  # Apr
    0.65,  # May
    0.66,  # Jun
    0.67,  # Jul
    0.66,  # Aug
    0.65,  # Sep
    0.66,  # Oct
    0.67,  # Nov
    0.66   # Dec
], dtype=float)


# Convert to sunshine fraction u = n/N ≈ 1 - cloud_fraction
sunshine_frac = 1.0 - CLOUD_FRAC_MONTH  # mean ≈ 0.335

# Your A–P-type mapping (keep as you had it)
def KT_from_sunshine(u): 
    return 0.18 + 0.62*float(u)

KT_month_raw  = np.array([KT_from_sunshine(u) for u in sunshine_frac])
KT_month_norm = KT_month_raw / KT_month_raw.mean()  # normalize so mean = 1
CLOUD_STRENGTH = 0.6  # 0..1: amplitude of cloud modulation
KT_month = 1.0 + CLOUD_STRENGTH*(KT_month_norm - 1.0)

# (Keep your existing day_to_month mapping)
day_to_month = np.repeat(np.arange(12), [31,28,31,30,31,30,31,31,30,31,30,31])[:365]


# --- Daisyworld core ---t
def clamp_cover(Aw, Ab, min_seed=0.01, max_total=0.99):
    Aw = max(min_seed, Aw); Ab = max(min_seed, Ab)
    s = Aw + Ab
    if s > max_total:
        Aw *= max_total/s; Ab *= max_total/s
    return Aw, Ab

def daily_step(Aw, Ab, L, I_eff, substeps=24):
    sub_dt = DT/substeps
    Tp=None
    for _ in range(substeps):
        x = max(0.0, 1.0 - Aw - Ab)
        Aplanet = Aw*ALB_W + Ab*ALB_B + x*ALB_G
        Tp = (L*I_eff*(1.0 - Aplanet)/SIGMA)**0.25
        Tw = ((R*L*I_eff/SIGMA*(Aplanet-ALB_W)) + Tp**4)**0.25
        Tb = ((R*L*I_eff/SIGMA*(Aplanet-ALB_B)) + Tp**4)**0.25
        betaw = 1.0 - 0.003265*(295.65 - Tw)**2 if T_MIN<=Tw<=T_MAX else 0.0
        betab = 1.0 - 0.003265*(295.65 - Tb)**2 if T_MIN<=Tb<=T_MAX else 0.0
        Aw += Aw*(x*betaw - GAMMA)*sub_dt
        Ab += Ab*(x*betab - GAMMA)*sub_dt
        Aw, Ab = clamp_cover(Aw, Ab)
    return Aw, Ab, Tp

def run_year(Aw, Ab, L):
    """Run one seasonal year with clouds; return end-state and annual-mean Tp."""
    Tp_year=[]
    for doy in range(1,366):
        m   = day_to_month[doy-1]
        KT  = KT_month[m]
        I_day = planetary_daily_mean(doy)
        I_eff = SCALE * I_day * KT
        Aw, Ab, Tp = daily_step(Aw, Ab, L, I_eff, substeps=24)
        Tp_year.append(Tp)
    return Aw, Ab, float(np.mean(Tp_year))

# --- Hysteresis sweep: carry state along L (up and down) ---
L_grid = np.arange(0.60, 1.80, 0.02)

# up-sweep
Aw_up=[]; Ab_up=[]; Tp_up=[]
Aw, Ab = 0.15, 0.08
for L in L_grid:
    # spin-up one year to equilibrate at this L
    Aw, Ab, _ = run_year(Aw, Ab, L)
    # measure one more year’s mean
    Aw, Ab, Tp_mean = run_year(Aw, Ab, L)
    Aw_up.append(Aw); Ab_up.append(Ab); Tp_up.append(Tp_mean)

# down-sweep (start from last up-state)
Aw_dn=[]; Ab_dn=[]; Tp_dn=[]
Aw, Ab = Aw_up[-1], Ab_up[-1]
for L in L_grid[::-1]:
    Aw, Ab, _ = run_year(Aw, Ab, L)
    Aw, Ab, Tp_mean = run_year(Aw, Ab, L)
    Aw_dn.append(Aw); Ab_dn.append(Ab); Tp_dn.append(Tp_mean)
Aw_dn = Aw_dn[::-1]; Ab_dn = Ab_dn[::-1]; Tp_dn = Tp_dn[::-1]

# bare-planet
def bare_T(L):
    I_mean = np.mean([planetary_daily_mean(d) for d in range(1,366)])
    return (L * SCALE * I_mean * (1.0 - ALB_G) / SIGMA)**0.25
Tp_bare = np.array([bare_T(L) for L in L_grid])



# fig, axs = plt.subplots(2,1, figsize=(7,10), sharex=True)

# axs[0].plot(L_grid, Tp_bare, color='gray', linestyle='--', label='Bare planet')
# axs[0].plot(L_grid, Tp_up, label='With daisies (up)')
# axs[0].plot(L_grid, Tp_dn, label='With daisies (down)')
# axs[0].set_ylabel('Mean Planet T (K)')
# axs[0].set_title('Equilibrium vs Luminosity (cloud-modulated, normalized)')
# axs[0].legend()

# axs[1].plot(L_grid, Aw_up, label='White (up)')
# axs[1].plot(L_grid, Ab_up, label='Black (up)')
# axs[1].plot(L_grid, Aw_dn, '--', label='White (down)')
# axs[1].plot(L_grid, Ab_dn, '--', label='Black (down)')
# axs[1].set_ylabel('Area fraction')
# axs[1].set_title('Daisy Cover vs Luminosity (hysteresis)')
# axs[1].legend(ncol=2)

plt.plot(range(1,13), (CLOUD_FRAC_MONTH), marker='o')
plt.axhline(y=0.66, color='red', linestyle='--', label='Mean cloud fraction 0.66')
plt.xlabel('Month')
plt.ylabel('Cloud Fraction $K_T$ (mean=0.65)')
plt.title('Monthly cloud modulation')

plt.tight_layout(); plt.show()
