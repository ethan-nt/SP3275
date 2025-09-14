'''
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

A_w_seq = []
A_b_seq = []
T_p_seq = []
time = []
gamma = 0.3
S_0 = 917
sigma = 5.670367e-8
timestep = 1
R = 0.2
albedo_w = 0.75
albedo_b = 0.25
albedo_g = 0.5
T_min = 278.15
T_max = 313.15
L_list=np.arange(0.4,1.7,0.01)
T_p2 = []
A_w2 = []
A_b2 = []

t = 0
A_w = 0.01
A_b = 0.01


for L in L_list:
    t = 0
    A_w = 0.01
    A_b = 0.01
    for i in range(50):
        x = 1 - A_w - A_b
        time.append(t)
        albedo_p = A_w * albedo_w + x * albedo_g + A_b * albedo_b
        T_p = (L * S_0 * (1 - albedo_p)/sigma)**0.25
        T_w = ((R * L * S_0/sigma * (albedo_p - albedo_w)) + T_p**4)**0.25
        T_b = ((R * L * S_0 / sigma * (albedo_p - albedo_b)) + T_p ** 4) ** 0.25
        T_p_seq.append(T_p)
        if T_b <= T_max and T_b >= T_min:
            beta_b = 1 - 0.003265 * (295.65 - T_b) ** 2
        else:
            beta_b = 0
        if T_w <= T_max and T_w >= T_min:
            beta_w = 1 - 0.003265 * (295.65 - T_w)**2
        else:
            beta_w = 0
        delta_A_w = A_w * (x * beta_w - gamma) * timestep
        delta_A_b = A_b * (x * beta_b - gamma) * timestep

        if A_b < 0.01:
            A_b = 0.01
        else:
            A_b = A_b + delta_A_b
        if A_w < 0.01:
            A_w = 0.01
        else:
            A_w = A_w + delta_A_w
        A_w_seq.append(A_w)
        A_b_seq.append(A_b)

        t += timestep
    T_p2.append(T_p_seq[-1])
    A_w2.append(A_w_seq[-1])
    A_b2.append(A_b_seq[-1])





fig, axs = plt.subplots(3, 1, figsize=(6, 8))

axs[0].plot(L_list, T_p2)
axs[0].set_title("Temperature against Luminosity")
axs[0].set_xlabel("Luminosity")
axs[0].set_ylabel("Temperature")

axs[1].plot(L_list, A_w2)
axs[1].set_title("White Daisy Area against Luminosity")
axs[1].set_xlabel("Luminosity")
axs[1].set_ylabel("White Daisy")

axs[2].plot(L_list, A_b2)
axs[2].set_title("Black Daisy Area against Luminosity")
axs[2].set_xlabel("Luminosity")
axs[2].set_ylabel("Black Daisy")

plt.tight_layout()
plt.show()

import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

A_w_seq = []
A_b_seq = []
T_p_seq = []
time = []
gamma = 0.3
S_0 = 917
sigma = 5.670367e-8
timestep = 0.1
R = 0.2
albedo_w = 0.75
albedo_b = 0.25
albedo_g = 0.5
T_min = 278.15
T_max = 313.15
L = 1.2
T_p2 = []
A_w2 = []
A_b2 = []

t = 0
A_w = 0.7
A_b = 0.3
albedo_p_seq = []
rate_w_seq = []
rate_b_seq = []
T_w_seq = []
T_b_seq = []

for i in range(500):
    x = 1 - A_w - A_b
    time.append(t)
    albedo_p = A_w * albedo_w + x * albedo_g + A_b * albedo_b
    albedo_p_seq.append(albedo_p)
    T_p = (L * S_0 * (1 - albedo_p)/sigma)**0.25
    T_w = ((R * L * S_0/sigma * (albedo_p - albedo_w)) + T_p**4)**0.25
    T_w_seq.append(T_w)
    T_b = ((R * L * S_0 / sigma * (albedo_p - albedo_b)) + T_p ** 4) ** 0.25
    T_b_seq.append(T_b)
    T_p_seq.append(T_p)
    if T_b <= T_max and T_b >= T_min:
        beta_b = 1 - 0.003265 * (295.65 - T_b) ** 2
    else:
        beta_b = 0
    if T_w <= T_max and T_w >= T_min:
        beta_w = 1 - 0.003265 * (295.65 - T_w)**2
    else:
        beta_w = 0
    delta_A_w = A_w * (x * beta_w - gamma) * timestep
    rate_w_seq.append(delta_A_w/timestep)
    delta_A_b = A_b * (x * beta_b - gamma) * timestep
    rate_b_seq.append(delta_A_b/timestep)
    A_b = A_b + delta_A_b
    A_w = A_w + delta_A_w

    A_w_seq.append(A_w)
    A_b_seq.append(A_b)
    t += timestep




fig, axs = plt.subplots(4, 1, figsize=(6, 8))

axs[0].plot(time, A_w_seq, label = 'WhitEe Daisy')
axs[0].plot(time, A_b_seq, label = 'Black Daisy')
axs[0].legend()
axs[0].set_title("Daisy Area against Time")
axs[0].set_xlabel("Time")
axs[0].set_ylabel("Daisy Area")

axs[1].plot(time, T_p_seq, label = 'Planet')
axs[1].plot(time, T_w_seq, label = 'White Daisy')
axs[1].plot(time, T_b_seq, label = 'Black Daisy')
axs[1].legend()
axs[1].set_title("Planet Temperature against Time")
axs[1].set_xlabel("time")
axs[1].set_ylabel("Planet Temperature")

axs[2].plot(time, albedo_p_seq)
axs[2].set_title("Planet Albedo against Time")
axs[2].set_xlabel("Time")
axs[2].set_ylabel("Planet Albedo")

axs[3].plot(time, rate_w_seq, label = 'White Daisy')
axs[3].plot(time, rate_b_seq, label = 'Black Daisy')
axs[3].legend()
axs[3].set_title("Rate of Area Changing against Time")
axs[3].set_xlabel("Time")
axs[3].set_ylabel("Rate of Changing")

plt.tight_layout()
plt.show()

class Planet:
  def __init__(self, position, mass, velocity):
    self.position = np.array(position,dtype=float) # input a tuple or array of any dimension
    self.mass = mass
    self.velocity = np.array(velocity,dtype=float)
    self.acc = np.zeros(len(position),dtype=float) # initial acc is 0 on all axes before being updated

  def update(self, delta, other_planet1, other_planet2):
    # calculate separation between home planet and planet 1
    separation_vector1 = other_planet1.position - self.position # returns a vector
    separation_magnitude1=np.linalg.norm(separation_vector1)

    separation_vector2 = other_planet2.position - self.position # returns a vector
    separation_magnitude2=np.linalg.norm(separation_vector2)

    # to calculate new acceleration of home planet
    self.acc = (G * other_planet1.mass)/(separation_magnitude1**3) * separation_vector1 + (G * other_planet2.mass)/(separation_magnitude2**3) * separation_vector2

    # calculate new velocity of home planet
    self.velocity += self.acc * delta # delta is change in time (in update)

    # calculate new position of home planet after a certain time step (delta)
    self.position += self.velocity * delta

    # to return a copy of the planet's original state
  def copy(self):
    p = Planet(self.position.copy(), self.mass, self.velocity.copy())
    p.acc = self.acc.copy()
    return p'
'''

'''
class planet():
    def __init__(self, mass, radius, position, velocity):
        self.radius = radius
        self.position = np.array(position,dtype=float) 
        self.mass = mass
        self.velocity = np.array(velocity,dtype=float)
        self.acc = np.zeros(len(position),dtype=float)

    def get_position(self):
        return self.position
    
    def get_radius(self):
        return self.radius
    
    def get_mass(self):
        return self.mass
    
    def get_velocity(self):
        return self.velocity
    
    def get_acc(self):
        return self.acc
    
    def update_acc(self, force):
        self.acc = force/ self.get_mass()

    def update_v(self, timestep):
        self.velocity += self.get_acc() * timestep

    def update_pos(self, timestep):
        self.position += self.get_velocity() * timestep
    
    def get_distance(self, planet2):
        return np.linalg.norm(self.get_position() - planet2.get_position())

    def update(self, planet2, timestep):
        grav_force = G * self.get_mass() * planet2.get_mass()/ np.linalg.norm(self.get_position() - planet2.get_position())**3 * (self.get_position() - planet2.get_position())
        self.update_acc(-grav_force)
        planet2.update_acc(grav_force)
        self.update_v(timestep)
        planet2.update_v(timestep)
        self.update_pos(timestep)
        planet2.update_pos(timestep)
'''

import numpy as np
import math
import matplotlib.pyplot as plt
K = 273.15
G = 6.674e-11
gamma = 0.3
S_0 = 917
sigma = 5.670367e-8
timestep = 1
R = 0.2
albedo_w = 0.75
albedo_b = 0.25
albedo_g = 0.5
T_min = 278.15
T_max = 313.15
solar_constant = 970


def solar_declination(time): 
    gamma = 2 * math.pi * ((time%365) - 1)/365
    return (0.006918
            - 0.399912*math.cos(gamma) + 0.070257*math.sin(gamma)
            - 0.006758*math.cos(2*gamma) + 0.000907*math.sin(2*gamma)
            - 0.002697*math.cos(3*gamma) + 0.001480*math.sin(3*gamma))
#http://dx.doi.org/10.1016/j.rser.2015.11.044



# def hour_angle(latitude, time):
#     return math.acos(-math.tan(latitude) * math.tan(solar_declination(time)))

# def hour_angle(latitude_rad, time):
#     # latitude_rad: TRUE geographic latitude in radians
#     delta = solar_declination(time)  # radians
#     x = -math.tan(latitude_rad) * math.tan(delta)

#     eps = 1e-12
#     if x >= 1 - eps:      # Sun below horizon all day (no sunrise)
#         return 0.0
#     elif x <= -1 + eps:   # Sun above horizon all day (midnight sun)
#         return math.pi
#     else:
#         x = max(-1.0, min(1.0, x))
#         return math.acos(x)


def eccentricity_corrector(time):
    return 1 + 0.033 * math.cos(2*math.pi * (time%365)/365)

def hour_angle_sunset(latitude, delta):
    # latitude and delta in radians
    cos_omega_s = -math.tan(latitude) * math.tan(delta)
    cos_omega_s = max(-1.0, min(1.0, cos_omega_s))
    return math.acos(cos_omega_s)

def daily_radiation(time, latitude):
    latitude = math.radians(latitude)
    E_0 = eccentricity_corrector(time)
    delta = solar_declination(time)
    omega_s = hour_angle_sunset(latitude, delta)
    # Daily mean insolation (zero at night)
    I_0 = (solar_constant / math.pi) * E_0 * (
        omega_s * math.sin(latitude) * math.sin(delta) +
        math.cos(latitude) * math.cos(delta) * math.sin(omega_s)
    )
    # Clamp to zero if sun never rises
    return max(I_0, 0)

def area_fraction(latitude):
   latitude = math.radians(latitude)
   return abs((math.sin(latitude + 1/180 * math.pi) - math.sin(latitude))/2)

T_p_seq = []
A_w_seq = []
A_b_seq = []
A_w = np.ones(180) * 0.8
A_b = np.ones(180) * 0.2
T_p = 20 + K
albedo_p = np.zeros(180)
areas = np.vectorize(area_fraction)(np.arange(0, 180))
for time in range(365*2):
    delta_A_b = np.zeros(180)
    delta_A_w = np.zeros(180)
    T = np.zeros(180)
    for lat in range(0, 180):
        x = np.ones(180) - A_w - A_b
        albedo_p[lat] = A_w[lat] * albedo_w + x[lat] * albedo_g + A_b[lat] * albedo_b
        I_0 = daily_radiation(time, lat)  # <-- use daily radiation
        area_i = areas[lat]
        T[lat] = (I_0 * (1 - albedo_p[lat])/sigma)**0.25
        # print(albedo_p, T_p_i, I_0)
        T_w = ((R * I_0 / sigma * (albedo_p[lat] - albedo_w)) + T[lat] ** 4) ** 0.25
        T_b = ((R * I_0 / sigma * (albedo_p[lat] - albedo_b)) + T[lat] ** 4) ** 0.25
        if T_b <= T_max and T_b >= T_min:
            beta_b = 1 - 0.003265 * (295.65 - T_b) ** 2
        else:
            beta_b = 0
        if T_w <= T_max and T_w >= T_min:
            beta_w = 1 - 0.003265 * (295.65 - T_w)**2
        else:
            beta_w = 0
        delta_A_w[lat] += area_i * A_w[lat] * (x[lat] * beta_w - gamma) * 1  # timestep = 1 day
        delta_A_b[lat] += area_i * A_b[lat] * (x[lat] * beta_b - gamma) * 1
    T_p_new = np.dot(T, areas)
    T_p = np.ones(180) * T_p_new
    A_w += delta_A_w
    A_w_seq.append(np.mean(A_w))
    A_b += delta_A_b
    A_b_seq.append(np.mean(A_b))
    T_p_seq.append(T_p_new)

plt.plot(A_w_seq, label = 'White')
plt.plot(A_b_seq, label = 'Black')
plt.show()
plt.plot(T_p_seq, label = 'Temperature')
plt.legend()
plt.show()