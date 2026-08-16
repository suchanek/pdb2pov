"""
Rotation and placement.

The degree-to-radian constants are the C's, and deliberately wrong: 57.29 and
57.298 are both truncations of 180/pi.  They are kept verbatim because the
scene geometry they produce is quoted downstream -- quiltwright's
``docs/pdb2pov.md`` pins crambin's camera distance at 40.075 and its enclosing
sphere at 18.759, both derived from these exact values.  Correcting them would
silently move every camera in every existing scene for a change of well under
a tenth of a degree.
"""

from __future__ import annotations

import math

from .structure import Structure

__all__ = [
    "DEG_PER_RAD_ROT",
    "DEG_PER_RAD_CAM",
    "make_rotmat",
    "matrix_times_vector",
    "rotate_structure",
]

DEG_PER_RAD_ROT = 57.29  #: used by :func:`make_rotmat`
DEG_PER_RAD_CAM = 57.298  #: used when placing the camera and the light

Matrix3 = tuple[
    tuple[float, float, float],
    tuple[float, float, float],
    tuple[float, float, float],
]


def make_rotmat(xdeg: float, ydeg: float, zdeg: float) -> Matrix3:
    """The C's rotation matrix, term for term."""
    g = xdeg / DEG_PER_RAD_ROT
    b = ydeg / DEG_PER_RAD_ROT
    t = zdeg / DEG_PER_RAD_ROT

    sing, cosg = math.sin(g), math.cos(g)
    sinb, cosb = math.sin(b), math.cos(b)
    sint, cost = math.sin(t), math.cos(t)

    return (
        (
            cost * cosb,
            sint * cosg + cost * sinb * sing,
            sint * sing - cost * sinb * cosg,
        ),
        (
            -sint * cosb,
            cost * cosg - sint * sinb * sing,
            cost * sing + sint * sinb * cosg,
        ),
        (
            sinb,
            -cosb * sing,
            cosb * cosg,
        ),
    )


def matrix_times_vector(
    m: Matrix3, x: float, y: float, z: float
) -> tuple[float, float, float]:
    return (
        m[0][0] * x + m[0][1] * y + m[0][2] * z,
        m[1][0] * x + m[1][1] * y + m[1][2] * z,
        m[2][0] * x + m[2][1] * y + m[2][2] * z,
    )


def rotate_structure(s: Structure, xdeg: float, ydeg: float, zdeg: float) -> None:
    """Apply absolute rotations about each axis, in place."""
    if xdeg == 0.0 and ydeg == 0.0 and zdeg == 0.0:
        return

    m = make_rotmat(xdeg, ydeg, zdeg)
    for a in s.atoms:
        a.x, a.y, a.z = matrix_times_vector(m, a.x, a.y, a.z)
