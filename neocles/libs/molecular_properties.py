import math
import inspect
from functools import wraps

import quantities as pq

pq.piconewton = pq.UnitQuantity('piconewton', (pq.newton * 1e-12).simplified,
                                symbol='pN')
pq.piconewtons = pq.piconewton
pq.kilojoule = pq.UnitQuantity('kilojoule', (pq.joule * 1e3).simplified,
                               symbol='kJ')
pq.viscosity.water = 1.00e-3 * pq.pascal * pq.second

T = 310.0 * pq.kelvin


class OutputUnitMisMatch(Exception):
    pass


def output_units(units):
    def wrapper(fn):
        @wraps(fn)
        def wrapped(*args, **kwargs):
            result = fn(*args, **kwargs)
            if result is not None and \
               result.dimensionality != units.dimensionality:
                raise OutputUnitMisMatch(fn.__name__, result, units)
            return result
        wrapped._args = inspect.getargspec(fn).args
        return wrapped
    return wrapper


# TODO: No rescale
@output_units(pq.kilogram)
def compute_mass(weight):
    return weight.rescale(pq.kilogram)


@output_units(pq.meter ** 3)
def compute_volume(mass, density):
    return mass / density


@output_units(pq.meter)
def compute_radius(volume):
    return (volume / (4.0/3.0 * math.pi)) ** (1.0/3.0)


@output_units(pq.newton * pq.second / pq.meter)
def compute_friction_coefficient(radius):
    value = 6.0 * math.pi * pq.viscosity.water * radius
    return value.rescale(pq.newton * pq.second / pq.meter)


@output_units(pq.piconewton * pq.second / pq.meter)
def compute_pico_friction_coefficient(friction_coefficient):
    return friction_coefficient.rescale(pq.piconewton * pq.second / pq.meter)


@output_units(pq.joule)
def compute_average_ke():
    value = 1.5 * pq.constants.Boltzmann_constant * T
    return value.rescale(pq.joule)


@output_units(pq.piconewton * pq.picometer)
def compute_average_pico_ke(average_ke):
    return average_ke.rescale(pq.piconewton * pq.picometer)


@output_units(pq.meter / pq.second)
def compute_average_velocity(mass):
    value = (3.0 * pq.constants.Boltzmann_constant * T / mass) ** 0.5
    return value.simplified


@output_units(pq.piconewton)
def compute_average_thermal_force(friction_coefficient, average_velocity):
    value = friction_coefficient * average_velocity
    return value.rescale(pq.piconewton)


@output_units(pq.piconewtons)
def compute_gravitational_force(mass):
    value = mass * pq.acceleration.gravity
    return value.rescale(pq.piconewton)


@output_units(pq.piconewton)
def compute_centrifugal_force(mass):
    value = mass * pq.acceleration.gravity * 100000.0
    return value.rescale(pq.piconewton)


@output_units(pq.meter ** 2 / pq.second)
def compute_linear_diffusion_coefficient(friction_coefficient):
    value = pq.constants.Boltzmann_constant * T / friction_coefficient
    return value.simplified


@output_units(pq.micrometer ** 2 / pq.second)
def compute_micro_linear_diffusion_coefficient(friction_coefficient):
    value = pq.constants.Boltzmann_constant * T / friction_coefficient
    return value.rescale(pq.micrometer ** 2 / pq.second)


@output_units(pq.picosecond)
def compute_translation_time_constant(mass, friction_coefficient):
    value = mass / friction_coefficient
    return value.rescale(pq.picosecond)


@output_units(pq.nanometer)
def compute_average_distance(average_velocity, translation_time_constant):
    value = (average_velocity * translation_time_constant)
    return value.rescale(pq.nanometer)


@output_units(pq.nanosecond)
def compute_first_passage_time(linear_diffusion_coefficient):
    value = (5.0 * pq.nanometer) ** 2 / (2.0 * linear_diffusion_coefficient)
    return value.rescale(pq.nanosecond)


@output_units(pq.nanometer)
def compute_diffusion_in_10_ps(linear_diffusion_coefficient):
    value = (10.0 * pq.picosecond * 2.0 * linear_diffusion_coefficient) ** 0.5
    return value.rescale(pq.nanometer)


@output_units(pq.nanometer)
def compute_diffusion_in_10_ns(linear_diffusion_coefficient):
    value = (10.0 * pq.nanosecond * 2.0 * linear_diffusion_coefficient) ** 0.5
    return value.rescale(pq.nanometer)


@output_units(pq.nanometer)
def compute_diffusion_in_10_us(linear_diffusion_coefficient):
    value = (10.0 * pq.microsecond * 2.0 * linear_diffusion_coefficient) ** 0.5
    return value.rescale(pq.nanometer)


@output_units(pq.nanometer)
def compute_diffusion_in_10_ms(linear_diffusion_coefficient):
    value = (10.0 * pq.millisecond * 2.0 * linear_diffusion_coefficient) ** 0.5
    return value.rescale(pq.nanometer)


@output_units(pq.meter / pq.second)
def compute_viscous_terminal_velocity(friction_coefficient,
                                      translation_time_constant):
    time = 1 * pq.nanosecond
    force = 1 * pq.piconewton
    exponent = -1 * (time / translation_time_constant).simplified
    value = (force / friction_coefficient) * (1.0 - math.exp(exponent))
    return value.simplified


@output_units(pq.nanometer)
def compute_viscous_distance(friction_coefficient,
                             translation_time_constant):
    time = 1 * pq.nanosecond
    force = 1 * pq.piconewton
    d1 = viscous_distance(time, force, friction_coefficient,
                          translation_time_constant)
    d0 = viscous_distance(0 * pq.picosecond, force, friction_coefficient,
                          translation_time_constant)
    return d1 - d0


def viscous_distance(time, force, friction_coefficient,
                     translation_time_constant):
    coeff = (force / friction_coefficient).simplified
    exp = (-1 * time / translation_time_constant).simplified
    value = coeff * (translation_time_constant * math.e ** exp + time)
    return value.rescale(pq.nanometer)


@output_units(pq.meter / pq.second)
def compute_ideal_terminal_velocity(mass):
    force = 1 * pq.piconewton
    time = 1 * pq.picosecond
    value = (force / mass) * time
    return value.simplified


@output_units(pq.nanometer)
def compute_ideal_distance(mass):
    force = 1 * pq.piconewton
    time = 1 * pq.picosecond
    value = 0.5 * (force / mass) * time ** 2
    return value.rescale(pq.nanometer)


@output_units(pq.joule)
def compute_thermal_energy():
    value = (1.0 / 2.0) * pq.constants.Boltzmann_constant * T
    return value.rescale(pq.joule)


@output_units(pq.piconewton * pq.nanometer)
def compute_molecule_thermal_energy(thermal_energy):
    return thermal_energy.rescale(pq.piconewton * pq.nanometer)


@output_units(pq.joule / pq.mole)
def compute_atp_per_mole():
    # CITE: http://en.wikipedia.org/wiki/ATP_hydrolysis
    return (-30.5 * pq.kilojoule / pq.mole).rescale(pq.joule / pq.mole)


@output_units(pq.piconewton * pq.nanometer)
def compute_atp_per_molecule(atp_per_mole):
    value = atp_per_mole * pq.mole
    return value.rescale(pq.piconewton * pq.nanometer)


@output_units(pq.dimensionless)
def compute_atp_kt_mulitples(atp_per_molecule, thermal_energy):
    return (atp_per_molecule / thermal_energy).simplified


@output_units(pq.joule / pq.mole)
def compute_h_bond_per_mol():
    # CITE: http://en.wikipedia.org/wiki/Hydrogen_bond
    return (5 * pq.kilojoule / pq.mole).rescale(pq.joule / pq.mole)


@output_units(pq.piconewton * pq.nanometer)
def compute_h_bond_per_molecule(h_bond_per_mol):
    value = h_bond_per_mol * pq.mole
    return value.rescale(pq.piconewton * pq.nanometer)


@output_units(pq.dimensionless)
def compute_h_bond_kt_mulitples(h_bond_per_molecule, thermal_energy):
    return (h_bond_per_molecule / thermal_energy).simplified


@output_units(pq.joule / pq.mole)
def compute_covalent_bond_energy():
    # CITE: http://www.science.uwaterloo.ca/~cchieh/cact/c120/bondel.html
    value = 348 * pq.kilojoule / pq.mole
    return value.rescale(pq.joule / pq.mole)


@output_units(pq.piconewton * pq.nanometer)
def compute_covalent_bond_per_molecule(covalent_bond_energy):
    value = covalent_bond_energy * pq.mole
    return value.rescale(pq.piconewton * pq.nanometer)


@output_units(pq.dimensionless)
def compute_covalent_bond_kt_mulitples(covalent_bond_per_molecule,
                                       thermal_energy):
    return (covalent_bond_per_molecule / thermal_energy).simplified


@output_units(pq.piconewton)
def compute_force_in_water():
    e = 80
    return electric_force(e, -1 * pq.elementary_charge, pq.elementary_charge,
                          pq.nanometer)


@output_units(pq.piconewton)
def compute_force_in_protein():
    e = 3
    return electric_force(e, -1 * pq.elementary_charge, pq.elementary_charge,
                          pq.nanometer)


@output_units(pq.piconewton)
def compute_force_in_vacuum():
    e = 1
    return electric_force(e, -1 * pq.elementary_charge, pq.elementary_charge,
                          pq.nanometer)


def electric_force(e, q_1, q_2, r):
    coeff = (1 / (4 * math.pi * pq.constants.electric_constant * e)).simplified
    ratio = (q_1 * q_2) / r ** 2
    return (coeff * ratio).rescale(pq.piconewton)
