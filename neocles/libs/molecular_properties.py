import math
import inspect
from functools import wraps

import quantities as pq

pq.piconewton = pq.UnitQuantity('piconewton', (pq.newton * 1e-12).simplified,
                                symbol='pN')
pq.piconewtons = pq.piconewton
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
    value = 0.5 * pq.constants.Boltzmann_constant * T
    return value.rescale(pq.joule)


@output_units(pq.piconewton * pq.picometer)
def compute_averge_pico_ke(average_ke):
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
def compute_viscous_distance(mass, friction_coefficient):
    time = 1 * pq.nanosecond
    force = 1 * pq.piconewton
    exp = (-1 * time * friction_coefficient / mass).simplified
    numerator = force * (mass * math.e ** exp + time * friction_coefficient)
    value = numerator / friction_coefficient ** 2
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
