# import numpy as np
import pandas as pd
import streamlit as st
import math
from scipy.constants import epsilon_0, speed_of_light, h, k
from arc import Lithium6, DynamicPolarizability

# Constants
mass = Lithium6.mass  # Mass of Lithium-6 in kg

# Streamlit app
st.title('Trap Parameters for Li6')

# Sidebar input fields
st.sidebar.header('Input Parameters')

power_slider = st.sidebar.slider('Laser Power, μW', min_value=1, max_value=1_000, value=80, step=1)
waist_slider = st.sidebar.slider('Beam Waist, μm', min_value=0.1, max_value=3.0, value=0.8, step=0.01)

# Wavelength
wavelength_input = st.sidebar.number_input('Wavelength, nm', value=808)
wavelength_input *= 1e-9  # Convert to meters

# Power
power_input = st.sidebar.number_input('Laser Power, μW', value=power_slider)
power_input *= 1e-6  # Convert to watts

# Waist
waist_input = st.sidebar.number_input('Beam Waist, μm', value=waist_slider)
waist_input *= 1e-6  # Convert to meters

# Cache the polarizability calculation to optimize performance
@st.cache_data
def get_alpha0(wavelength, n, l, j):
    atom = Lithium6(preferQuantumDefects=False, cpp_numerov=True)
    calc = DynamicPolarizability(atom, n, l, j)
    calc.defineBasis(nMin=atom.groundStateN, nMax=n+15)
    alpha0, *_ = calc.getPolarizability(wavelength, units='SI', accountForStateLifetime=True)
    return alpha0

# Calculate U_0 (Trap depth)
def get_U0(wavelength, power, waist):
    alpha0 = get_alpha0(wavelength, n=2, l=0, j=1/2)
    return alpha0 * power / (epsilon_0 * speed_of_light * math.pi * waist**2)

# Calculate deltaE (Excitation energy)
def get_deltaE(wavelength, power, waist):
    alpha0_2S = get_alpha0(wavelength, n=2, l=0, j=1/2)
    alpha0_2P = get_alpha0(wavelength, n=2, l=1, j=3/2)
    return power / (epsilon_0 * speed_of_light * math.pi * waist**2) * (alpha0_2S + abs(alpha0_2P))

# Calculate radial trap frequency (omega_r)
def get_omega_r(U0, waist):
    return math.sqrt(4 * U0 / (mass * waist**2))

# Calculate axial trap frequency (omega_z)
def get_omega_z(U0, wavelength, waist):
    return math.sqrt(2 * U0 * wavelength**2 / (mass * math.pi**2 * waist**4))

def get_omega_ratio(wavelength, power, waist):
    U0 = get_U0(wavelength, power, waist) * h
    return get_omega_r(U0, waist) / get_omega_z(U0, wavelength, waist)

# Calculate and display results
U0 = get_U0(wavelength_input, power_input, waist_input) * h  # U0 in joules
deltaE = get_deltaE(wavelength_input, power_input, waist_input)

# Radial and axial trap frequencies
omega_r = get_omega_r(U0, waist_input)
omega_z = get_omega_z(U0, wavelength_input, waist_input)

# 
on = st.sidebar.toggle("Switch to temperature")

# Create a dataframe for displaying results
if not on:
    data = pd.DataFrame({
        "Parameter": ["Trap depth", "D2 detuning", "Trap frequency: radial", "Trap frequency: axial", "Trap frequency: ratio"],
        "Symbol": ["U₀", "ΔE", "ω_r / 2π", "ω_z / 2π", "ω_r / ω_z"],
        "Value": [U0 / h * 1e-6, deltaE * 1e-6, omega_r / (2 * math.pi) * 1e-3, omega_z / (2 * math.pi) * 1e-3, omega_r / omega_z],
        "Units": ["MHz", "MHz", "kHz", "kHz", ""]
    })
else:
    data = pd.DataFrame({
        "Parameter": ["Trap depth", "D2 detuning", "Trap frequency: radial", "Trap frequency: axial", "Trap frequency: ratio"],
        "Symbol": ["U₀", "ΔE", "ω_r / 2π", "ω_z / 2π", "ω_r / ω_z"],
        "Value": [U0 / k * 1e6, deltaE * 1e-6, omega_r / (2 * math.pi) * 1e-3, omega_z / (2 * math.pi) * 1e-3, omega_r / omega_z],
        "Units": ["μK", "MHz", "kHz", "kHz", ""]
    })

st.write(data)