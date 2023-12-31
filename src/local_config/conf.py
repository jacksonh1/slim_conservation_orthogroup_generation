# from enum import Enum
# class MyEnum(Enum):
#     """Docstring for MyEnum."""
#     FIRST_ENUM = "some_value"
#     SECOND_ENUM = "some_other_value"
    
import copy

from attrs import define, field, validators

import local_env_variables.env_variables as env


@define
class FilterConf:
    """sequence filtering parameters"""    
    min_fraction_short_than_query: float = field(
        default=0.5,
        validator=validators.and_(validators.le(1), validators.ge(0))
    )

@define
class OGSelectConf:
    OG_selection_method: str = field(
        default="level_name",
        validator=validators.in_(["most_species", "level_name"])
    )
    OG_level_name: str = field(default="Vertebrata")

@define
class LDOSelectConf:
    LDO_selection_method: str = field(
        default="alfpy_google_distance",
        validator=validators.in_(
            ["msa_by_organism", "alfpy_google_distance", "pairwise", "msa"]
        )
    )

@define
class AlignConf:
    align: bool = field(default=False, converter=bool)
    n_align_threads: int = field(default=8, converter=int)
    mafft_exe: str = field(default=env.MAFFT_EXECUTABLE)
    mafft_additional_args: str = field(default=env.MAFFT_ADDITIONAL_ARGUMENTS)

@define
class PipelineParams:
    filter_params: FilterConf = field(default=FilterConf())
    og_select_params: OGSelectConf = field(default=OGSelectConf())
    ldo_select_params: LDOSelectConf = field(default=LDOSelectConf())
    align_params: AlignConf = field(default=AlignConf())
    cd_hit_exe: str = field(default=env.CD_HIT_EXECUTABLE)
    cd_hit_additional_args: str = field(default=env.CD_HIT_ADDITIONAL_ARGUMENTS)
    main_output_folder: str = field(default="./odb_group_construction_output")
    write_files: bool = field(default=True)

    @classmethod
    def from_dict(cls, d):
        d = copy.deepcopy(d)
        return cls(
            filter_params=FilterConf(**d.pop("filter_params", {})),
            og_select_params=OGSelectConf(**d.pop("og_select_params", {})),
            ldo_select_params=LDOSelectConf(**d.pop("ldo_select_params", {})),
            align_params=AlignConf(**d.pop("align_params", {})),
            **d,
        )


# ==============================================================================
# // test
# ==============================================================================

# DEFAULT_PARAM_DICT = {
#     "filter_params": {
#         "min_fraction_short_than_query": 0.5,
#     },
#     "og_select_params": {
#         "OG_selection_method": "level_name",
#         "OG_level_name": "Eukaryota",
#     },
#     "ldo_select_params": {
#         "LDO_selection_method": "alfpy_google_distance",
#     },
#     "align_params": {
#         "align": False,
#         "n_align_threads": 8,
#     },
#     "write_files": True,
#     "main_output_folder": "./orthoDB_analysis",
# }
# DEFAULT_PARAM_DICT.pop

# config = PipelineParams.from_dict(DEFAULT_PARAM_DICT)
# print(config)
# print(config.filter_params.min_fraction_short_than_query)
# print(DEFAULT_PARAM_DICT)







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
