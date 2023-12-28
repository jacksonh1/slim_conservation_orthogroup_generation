
import copy

from attrs import define, field, validators


@define
class filter_conf:
    min_fraction_short_than_query: float = field(default=0.5)

@define
class OG_select_conf:
    OG_selection_method: str = field(
        default="level name",
        validator=validators.in_(["most_species", "level name"])
    )
    OG_level_name: str = field(default="Eukaryota")


@define
class LDO_select_conf:
    LDO_selection_method: str = field(
        default="alfpy_google_distance",
        validator=validators.in_(
            ["msa_by_organism", "alfpy_google_distance", "pairwise", "msa"]
        )
    )


@define
class align_conf:
    align: bool = field(default=False, converter=bool)
    n_align_threads: int = field(default=32, converter=int)


@define
class pipeline_params:
    filter_params: filter_conf
    og_select_params: OG_select_conf
    ldo_select_params: LDO_select_conf
    align_params: align_conf
    main_output_folder: str = field(default="./orthoDB_analysis")
    write_files: bool = field(default=True)

    @classmethod
    def from_dict(cls, d):
        d = copy.deepcopy(d)
        return cls(
            filter_params=filter_conf(**d.pop("filter_params", {})),
            og_select_params=OG_select_conf(**d.pop("og_select_params", {})),
            ldo_select_params=LDO_select_conf(**d.pop("ldo_select_params", {})),
            align_params=align_conf(**d.pop("align_params", {})),
            **d,
        )


# ==============================================================================
# // test
# ==============================================================================

DEFAULT_PARAM_DICT = {
    "filter_params": {
        "min_fraction_short_than_query": 0.5,
    },
    "og_select_params": {
        "OG_selection_method": "level name",
        "OG_level_name": "Eukaryota",
    },
    # "ldo_select_params": {
        # "LDO_selection_method": "alfpy_google_distance",
    # },
    "align_params": {
        "align": False,
        "n_align_threads": 32,
    },
    "write_files": True,
    "main_output_folder": "./orthoDB_analysis",
}
# DEFAULT_PARAM_DICT.pop

config = pipeline_params.from_dict(DEFAULT_PARAM_DICT)
print(config)
print(config.filter_params.min_fraction_short_than_query)
print(DEFAULT_PARAM_DICT)







# @define
# class pipeline_conf:
#     main_output_folder: str = field(default="./orthoDB_analysis")
#     OG_selection_method: str = field(default="level name")
#     min_fraction_short_than_query: float = field(default=0.5)
#     OG_level_name: str = field(default="Eukaryota")
#     LDO_selection_method: str = field(default="alfpy_google_distance")
#     align: bool = field(default=False)
#     n_align_threads: int = field(default=32)
#     write_files: bool = field(default=True)
