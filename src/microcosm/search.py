"""Objective-value helpers for a community-FBA result.

`target_value` reads the net secretion (produce mode) or net uptake (consume mode -- for greenhouse-gas or
pollutant capture, where the aim is to draw a species down) of a named target from a
`community.community_transform` result. These are the small, pure adapters that turn a transform into the
scalar an analysis optimises or compares.
"""
from __future__ import annotations


def target_value(result, target, consume=False):
    """The objective magnitude for `target`: net SECRETION (produce mode) or net UPTAKE (consume mode -- for
    greenhouse-gas / pollutant capture, where the aim is to draw the target DOWN). 0 if the sign is wrong."""
    v = result.get("transform", {}).get(target, 0.0)
    return max(0.0, -v) if consume else max(0.0, v)


def target_production(result, target):                 # produce-mode convenience wrapper
    return target_value(result, target, consume=False)
