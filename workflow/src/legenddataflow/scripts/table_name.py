from legendmeta import LegendMetadata


def get_table_name(metadata_path, timestamp, datatype, detector):
    meta = LegendMetadata(path=metadata_path)
    channel_dict = meta.channelmap(timestamp, system=datatype)
    return f"ch{channel_dict[detector].daq.rawid:07}"
