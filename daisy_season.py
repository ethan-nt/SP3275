'''
import numpy as np
import math
import matplotlib.pyplot as plt

# -----------------------------
# Constants and parameters
# -----------------------------
GAMMA = 0.3                # death rate
S0 = 917.0                 # stellar constant used in your model (W/m^2) — keep as in your code
SIGMA = 5.670367e-8        # Stefan-Boltzmann constant (W m^-2 K^-4)
R = 0.2                    # local heating/cooling coupling
ALB_W, ALB_B, ALB_G = 0.75, 0.25, 0.50
T_MIN, T_MAX = 278.15, 313.15
DT = 1.0                   # time step for the daisy dynamics (in "model" units)

# growth-rate shape (centered near 295.65 K, as in your code)
def growth(T):
    if (T < T_MIN) or (T > T_MAX):
        return 0.0
    return 1.0 - 0.003265*(295.65 - T)**2

def clamp_cover(Aw, Ab, min_seed=0.01, max_total=0.99):
    # keep minimum seeds and cap total coverage
    Aw = max(min_seed, Aw)
    Ab = max(min_seed, Ab)
    s = Aw + Ab
    if s > max_total:
        Aw *= max_total/s
        Ab *= max_total/s
    return Aw, Ab

# One Euler step of Daisyworld at global-mean insolation I (W/m^2) with luminosity L
def step_once(L, Aw, Ab, I):
    x = max(0.0, 1.0 - Aw - Ab)  # bare ground
    # planetary albedo
    Aplanet = Aw*ALB_W + Ab*ALB_B + x*ALB_G
    # planetary temperature
    Tp = (L * I * (1.0 - Aplanet) / SIGMA)**0.25

    # local daisy temperatures
    Tw = ((R * L * I / SIGMA * (Aplanet - ALB_W)) + Tp**4)**0.25
    Tb = ((R * L * I / SIGMA * (Aplanet - ALB_B)) + Tp**4)**0.25

    # growth rates
    betaw = growth(Tw)
    betab = growth(Tb)

    # area changes
    dAw = Aw*(x*betaw - GAMMA)*DT
    dAb = Ab*(x*betab - GAMMA)*DT

    Aw_new = Aw + dAw
    Ab_new = Ab + dAb
    Aw_new, Ab_new = clamp_cover(Aw_new, Ab_new)
    return Aw_new, Ab_new, Tp, Aplanet, Tw, Tb

def converge_at_L(L, Aw0, Ab0, I, tol=1e-6, max_iter=5000):
    Aw, Ab = Aw0, Ab0
    Tp_prev = None
    for _ in range(max_iter):
        Aw_new, Ab_new, Tp, _, _, _ = step_once(L, Aw, Ab, I)
        if Tp_prev is not None and abs(Tp - Tp_prev) < tol and abs(Aw_new - Aw) < tol and abs(Ab_new - Ab) < tol:
            return Aw_new, Ab_new, Tp
        Aw, Ab, Tp_prev = Aw_new, Ab_new, Tp
    return Aw, Ab, Tp  # return last state if not fully converged

# -----------------------------
# Part A: Temperature vs Luminosity (equilibria)
# -----------------------------
L_list = np.arange(0.4, 1.7+1e-9, 0.01)

# For this sweep we use a fixed global-mean insolation "I" (your S0) to match your first block
I_global = S0

Aw_up, Ab_up, Tp_up = [], [], []
Aw, Ab = 0.01, 0.01  # seed once
for L in L_list:
    Aw, Ab, Tp = converge_at_L(L, Aw, Ab, I_global)
    Aw_up.append(Aw); Ab_up.append(Ab); Tp_up.append(Tp)

# (Optional) down-sweep to see hysteresis — start from last equilibrium
Aw_dn, Ab_dn, Tp_dn = [], [], []
Aw, Ab = Aw_up[-1], Ab_up[-1]
for L in L_list[::-1]:
    Aw, Ab, Tp = converge_at_L(L, Aw, Ab, I_global)
    Aw_dn.append(Aw); Ab_dn.append(Ab); Tp_dn.append(Tp)
Aw_dn, Ab_dn, Tp_dn = Aw_dn[::-1], Ab_dn[::-1], Tp_dn[::-1]

fig1, axs = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
axs[0].plot(L_list, Tp_up, label='Up-sweep')
axs[0].plot(L_list, Tp_dn, '--', label='Down-sweep')
axs[0].set_ylabel('Planet T (K)')
axs[0].set_title('Planet Temperature vs Luminosity')
axs[0].legend()

axs[1].plot(L_list, Aw_up, label='White (up)')
axs[1].plot(L_list, Ab_up, label='Black (up)')
axs[1].set_ylabel('Area fraction')
axs[1].set_title('Daisy Cover vs Luminosity (up-sweep)')
axs[1].legend()

axs[2].plot(L_list, Aw_dn, label='White (down)')
axs[2].plot(L_list, Ab_dn, label='Black (down)')
axs[2].set_xlabel('Luminosity (L)')
axs[2].set_ylabel('Area fraction')
axs[2].set_title('Daisy Cover vs Luminosity (down-sweep)')
axs[2].legend()

plt.tight_layout()
plt.show()

# -----------------------------
# Part B: Temperature vs Time within a Year (seasonal)
# -----------------------------
# We'll compute daily-mean TOA insolation as a function of latitude and day-of-year,
# then area-average over latitude bands to get a planetary-mean insolation for the model.

S0_TOA = 1361.0  # physical solar constant at TOA (W/m^2) for seasonal geometry
def deg2rad(d): return d * math.pi / 180.0

def declination(day):
    # Simple approximation (radians)
    return deg2rad(23.44) * math.sin(2*math.pi*(day-81)/365.0)

def eccentricity_correction(day):
    return 1.0 + 0.033 * math.cos(2*math.pi*day/365.0)

def daily_mean_toa(lat_deg, day):
    phi = deg2rad(lat_deg)
    delta = declination(day)
    E0 = eccentricity_correction(day)
    # sunset hour angle (clamped for polar day/night)
    arg = -math.tan(phi)*math.tan(delta)
    arg = max(-1.0, min(1.0, arg))
    omega_s = math.acos(arg)
    return (24.0/math.pi) * S0_TOA * E0 * (
        omega_s * math.sin(phi) * math.sin(delta) +
        math.cos(phi) * math.cos(delta) * math.sin(omega_s)
    )

# Planetary mean daily insolation by latitudinal area-weighted average
lat_edges = np.linspace(-90, 90, 181)  # 1-degree bands
lat_centers = 0.5*(lat_edges[:-1] + lat_edges[1:])
# fractional area weight for each band on a sphere is proportional to cos(latitude)
weights = np.cos(np.deg2rad(lat_centers))
weights /= weights.sum()

def planetary_daily_mean(day):
    vals = np.array([daily_mean_toa(lat, day) for lat in lat_centers])
    return np.sum(vals * weights)

# Run one model year with seasonal forcing, at some chosen luminosity L_seasonal
L_seasonal = 1.20
days = np.arange(1, 366)  # 1..365

Aw_series, Ab_series = [], []
Tp_series = []

Aw, Ab = 0.5, 0.2  # start with some cover (you can change)
for d in days:
    # global-mean insolation for this day (W/m^2)
    I_day = planetary_daily_mean(d)

    # take a few daisy-dynamics substeps per day to relax on that day's forcing
    substeps = 8
    sub_dt = DT / substeps
    for _ in range(substeps):
        # temporarily set smaller DT for stability
        x = max(0.0, 1.0 - Aw - Ab)
        Aplanet = Aw*ALB_W + Ab*ALB_B + x*ALB_G
        Tp = (L_seasonal * I_day * (1.0 - Aplanet) / SIGMA)**0.25
        Tw = ((R * L_seasonal * I_day / SIGMA * (Aplanet - ALB_W)) + Tp**4)**0.25
        Tb = ((R * L_seasonal * I_day / SIGMA * (Aplanet - ALB_B)) + Tp**4)**0.25
        betaw, betab = growth(Tw), growth(Tb)
        dAw = Aw*(x*betaw - GAMMA)*sub_dt
        dAb = Ab*(x*betab - GAMMA)*sub_dt
        Aw, Ab = Aw + dAw, Ab + dAb
        Aw, Ab = clamp_cover(Aw, Ab)

    # record end-of-day values and the consistent Tp for that day
    x = max(0.0, 1.0 - Aw - Ab)
    Aplanet = Aw*ALB_W + Ab*ALB_B + x*ALB_G
    Tp = (L_seasonal * I_day * (1.0 - Aplanet) / SIGMA)**0.25

    Aw_series.append(Aw); Ab_series.append(Ab); Tp_series.append(Tp)

fig2, axs2 = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
axs2[0].plot(days, Tp_series)
axs2[0].set_ylabel('Planet T (K)')
axs2[0].set_title(f'Planet Temperature vs Time (L = {L_seasonal:.2f})')

axs2[1].plot(days, Aw_series, label='White')
axs2[1].plot(days, Ab_series, label='Black')
axs2[1].set_ylabel('Area fraction')
axs2[1].set_title('Daisy Cover vs Time')
axs2[1].legend()

# (optional) show the planetary-mean insolation used each day
I_series = [planetary_daily_mean(d) for d in days]
axs2[2].plot(days, I_series)
axs2[2].set_xlabel('Day of year')
axs2[2].set_ylabel('TOA daily-mean (W/m²)')
axs2[2].set_title('Daily-mean TOA Insolation (planetary mean)')

plt.tight_layout()
plt.show()

import numpy as np
import math
import matplotlib.pyplot as plt

# -----------------------------
# Constants
# -----------------------------
GAMMA = 0.3                 # death rate
S0 = 917.0                  # stellar constant used in Daisyworld model (W/m^2)
SIGMA = 5.670367e-8         # Stefan–Boltzmann constant (W m^-2 K^-4)
R = 0.2                     # local heating/cooling coupling
ALB_W, ALB_B, ALB_G = 0.75, 0.25, 0.50
T_MIN, T_MAX = 278.15, 313.15
DT = 1.0                    # timestep for daisy dynamics

# Solar constant at Earth orbit for insolation function
S0_TOA = 1361.0             # W/m^2

# -----------------------------
# Helper functions
# -----------------------------
def clamp_cover(Aw, Ab, min_seed=0.01, max_total=0.99):
    """Keep min seed and cap total daisy cover"""
    Aw = max(min_seed, Aw)
    Ab = max(min_seed, Ab)
    s = Aw + Ab
    if s > max_total:
        Aw *= max_total/s
        Ab *= max_total/s
    return Aw, Ab

def deg2rad(d): 
    return d * math.pi / 180.0

def declination(day):
    """Solar declination angle δ in radians"""
    return deg2rad(23.44) * math.sin(2*math.pi*(day-81)/365.0)

def eccentricity_correction(day):
    """Earth–Sun distance correction"""
    return 1.0 + 0.033 * math.cos(2*math.pi*day/365.0)

def daily_mean_toa(lat_deg, day):
    """
    Daily-mean TOA irradiance (W/m^2) at latitude 'lat_deg'.
    Correct expression uses (S0/pi), not (24/pi).
    """
    phi = deg2rad(lat_deg)
    delta = declination(day)
    E0 = eccentricity_correction(day)

    # Sunset hour angle; clamp for polar day/night
    arg = -math.tan(phi) * math.tan(delta)
    arg = max(-1.0, min(1.0, arg))
    omega_s = math.acos(arg)

    return (S0_TOA / math.pi) * E0 * (
        omega_s * math.sin(phi) * math.sin(delta) +
        math.cos(phi) * math.cos(delta) * math.sin(omega_s)
    )

# Planetary mean daily insolation (area-weighted average)
lat_edges   = np.linspace(-90, 90, 181)
lat_centers = 0.5 * (lat_edges[:-1] + lat_edges[1:])
weights     = np.cos(np.deg2rad(lat_centers))
weights    /= weights.sum()

def planetary_daily_mean(day):
    vals = np.array([daily_mean_toa(lat, day) for lat in lat_centers])
    return np.sum(vals * weights)  # ≈ S0_TOA*E0/4

# -----------------------------
# Seasonal simulation
# -----------------------------
L_seasonal = 1.20
days = np.arange(1, 366)

Aw_series, Ab_series, Tp_series = [], [], []
Aw, Ab = 0.5, 0.2   # initial cover

for d in days:
    I_day = planetary_daily_mean(d)

    # do substeps for stability
    substeps = 8
    sub_dt = DT / substeps
    for _ in range(substeps):
        x = max(0.0, 1.0 - Aw - Ab)
        Aplanet = Aw*ALB_W + Ab*ALB_B + x*ALB_G
        Tp = (L_seasonal * I_day * (1.0 - Aplanet) / SIGMA)**0.25
        Tw = ((R * L_seasonal * I_day / SIGMA * (Aplanet - ALB_W)) + Tp**4)**0.25
        Tb = ((R * L_seasonal * I_day / SIGMA * (Aplanet - ALB_B)) + Tp**4)**0.25

        # growth rates
        betaw = 1.0 - 0.003265*(295.65 - Tw)**2 if T_MIN <= Tw <= T_MAX else 0.0
        betab = 1.0 - 0.003265*(295.65 - Tb)**2 if T_MIN <= Tb <= T_MAX else 0.0

        # update daisy fractions
        dAw = Aw*(x*betaw - GAMMA)*sub_dt
        dAb = Ab*(x*betab - GAMMA)*sub_dt
        Aw, Ab = Aw + dAw, Ab + dAb
        Aw, Ab = clamp_cover(Aw, Ab)

    # record daily values
    x = max(0.0, 1.0 - Aw - Ab)
    Aplanet = Aw*ALB_W + Ab*ALB_B + x*ALB_G
    Tp = (L_seasonal * I_day * (1.0 - Aplanet) / SIGMA)**0.25
    Aw_series.append(Aw)
    Ab_series.append(Ab)
    Tp_series.append(Tp)

# also store the insolation
I_series = [planetary_daily_mean(d) for d in days]

# -----------------------------
# Plot
# -----------------------------
fig, axs = plt.subplots(3, 1, figsize=(7, 9), sharex=True)

axs[0].plot(days, Tp_series)
axs[0].set_ylabel('Planet T (K)')
axs[0].set_title(f'Planet Temperature vs Time (L = {L_seasonal:.2f})')

axs[1].plot(days, Aw_series, label='White')
axs[1].plot(days, Ab_series, label='Black')
axs[1].set_ylabel('Area fraction')
axs[1].set_title('Daisy Cover vs Time')
axs[1].legend()

axs[2].plot(days, I_series)
axs[2].set_xlabel('Day of year')
axs[2].set_ylabel('TOA daily-mean (W/m²)')
axs[2].set_title('Daily-mean TOA Insolation (planetary mean)')

plt.tight_layout()
plt.show()
'''
# === Daisyworld: seasonal run (stand-alone) ===
import numpy as np
import math
import matplotlib.pyplot as plt

# -----------------------------
# Model constants (you can tweak)
# -----------------------------
SIGMA   = 5.670367e-8     # Stefan–Boltzmann (W m^-2 K^-4)
S0_TOA  = 1361.0          # Solar constant at TOA (W/m^2) for geometry
S0_MODEL= 917.0           # Daisyworld's internal stellar constant (your earlier code)
SCALE   = S0_MODEL / (S0_TOA / 4.0)  # map planetary-mean TOA → Daisyworld units (~2.697)

# Ecology / thermodynamics
ALB_W, ALB_B, ALB_G = 0.75, 0.25, 0.50   # albedos: white, black, bare ground
GAMMA  = 0.25          # death rate (slightly friendlier than 0.30)
R      = 0.25          # local temperature coupling (was 0.20)
T_MIN, T_MAX = 278.15, 313.15  # viable growth window (K)
DT     = 1.0           # model day (we’ll substep this)

# Luminosity used for the seasonal run (in Daisyworld units)
L_seasonal = 1.10

# -----------------------------
# Helper functions
# -----------------------------
def clamp_cover(Aw, Ab, min_seed=0.01, max_total=0.99):
    """Keep small seeds alive and cap total cover."""
    Aw = max(min_seed, Aw)
    Ab = max(min_seed, Ab)
    s = Aw + Ab
    if s > max_total:
        Aw *= max_total / s
        Ab *= max_total / s
    return Aw, Ab

def deg2rad(d): 
    return d * math.pi / 180.0

def declination(day):
    """Solar declination δ (rad), simple approximation."""
    return deg2rad(23.44) * math.sin(2 * math.pi * (day - 81) / 365.0)

def eccentricity_correction(day):
    """Seasonal 1/r^2 correction (dimensionless)."""
    return 1.0 + 0.033 * math.cos(2 * math.pi * day / 365.0)

def daily_mean_toa(lat_deg, day):
    """
    Daily-mean TOA irradiance (W/m^2) at latitude (Duffie & Beckman / Iqbal form).
    Note the (S0/pi) factor gives mean POWER (not daily total).
    """
    phi   = deg2rad(lat_deg)
    delta = declination(day)
    E0    = eccentricity_correction(day)
    # Sunset hour angle (clamped for polar conditions)
    arg = -math.tan(phi) * math.tan(delta)
    arg = max(-1.0, min(1.0, arg))
    omega_s = math.acos(arg)
    return (S0_TOA / math.pi) * E0 * (
        omega_s * math.sin(phi) * math.sin(delta) +
        math.cos(phi) * math.cos(delta) * math.sin(omega_s)
    )

# Area-weighted planetary mean (weights ∝ cos φ)
lat_edges   = np.linspace(-90, 90, 181)
lat_centers = 0.5 * (lat_edges[:-1] + lat_edges[1:])
weights     = np.cos(np.deg2rad(lat_centers))
weights    /= weights.sum()

def planetary_daily_mean(day):
    vals = np.array([daily_mean_toa(lat, day) for lat in lat_centers])
    return np.sum(vals * weights)   # ≈ S0_TOA*E0/4

# -----------------------------
# Seasonal simulation
# -----------------------------
days = np.arange(1, 366)

# Initial daisy cover (choose something non-tiny to avoid early extinction)
Aw, Ab = 0.1, 0.10

Aw_series, Ab_series, Tp_series = [], [], []
substeps = 24               # integrate daisies with smaller dt for stability
sub_dt  = DT / substeps

for d in days:
    I_day  = planetary_daily_mean(d)     # ~330–350 W/m^2
    I_eff  = SCALE * I_day               # map to Daisyworld's 917-based units

    for _ in range(substeps):
        x = max(0.0, 1.0 - Aw - Ab)                         # bare ground fraction
        Aplanet = Aw*ALB_W + Ab*ALB_B + x*ALB_G             # planetary albedo
        Tp = (L_seasonal * I_eff * (1.0 - Aplanet) / SIGMA)**0.25

        # Local daisy temperatures (standard Daisyworld parameterization)
        Tw = ((R * L_seasonal * I_eff / SIGMA * (Aplanet - ALB_W)) + Tp**4)**0.25
        Tb = ((R * L_seasonal * I_eff / SIGMA * (Aplanet - ALB_B)) + Tp**4)**0.25

        # Quadratic growth law centered near 295.65 K
        betaw = 1.0 - 0.003265 * (295.65 - Tw)**2 if (T_MIN <= Tw <= T_MAX) else 0.0
        betab = 1.0 - 0.003265 * (295.65 - Tb)**2 if (T_MIN <= Tb <= T_MAX) else 0.0

        # Area fraction tendencies
        dAw = Aw * (x * betaw - GAMMA) * sub_dt
        dAb = Ab * (x * betab - GAMMA) * sub_dt

        # Update and keep within bounds
        Aw += dAw; Ab += dAb
        Aw, Ab = clamp_cover(Aw, Ab)

    # Record end-of-day values
    x = max(0.0, 1.0 - Aw - Ab)
    Aplanet = Aw*ALB_W + Ab*ALB_B + x*ALB_G
    Tp = (L_seasonal * I_eff * (1.0 - Aplanet) / SIGMA)**0.25

    Aw_series.append(Aw)
    Ab_series.append(Ab)
    Tp_series.append(Tp)

# For plot 3: the actual forcing used each day (after scaling)
I_series = [SCALE * planetary_daily_mean(d) for d in days]

# -----------------------------
# Plots
# -----------------------------
fig, axs = plt.subplots(3, 1, figsize=(7, 9), sharex=True)

axs[0].plot(days, Tp_series)
axs[0].set_ylabel('Planet T (K)')
axs[0].set_title(f'Planet Temperature vs Time (L = {L_seasonal:.2f})')

axs[1].plot(days, Aw_series, label='White')
axs[1].plot(days, Ab_series, label='Black')
axs[1].set_ylabel('Area fraction')
axs[1].set_title('Daisy Cover vs Time')
axs[1].legend()

axs[2].plot(days, I_series)
axs[2].set_xlabel('Day of year')
axs[2].set_ylabel('Forcing $I_\\mathrm{{eff}}$ (W m$^{{-2}}$)')
axs[2].set_title('Planetary-mean Forcing (rescaled to model units)')

plt.tight_layout()
plt.show()
