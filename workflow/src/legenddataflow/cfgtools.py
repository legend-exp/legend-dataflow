from typing import Mapping


def get_channel_config(
    mapping: Mapping, channel: str, default_key: str = "__default__"
):
    """Get channel key from mapping with default.

    Returns the value at key `channel`, if existing, otherwise return value at
    `default_key`.
    """
    return mapping.get(channel, mapping[default_key])
