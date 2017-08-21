#include "gPMC_ggo.h"

template < class TArgsInfo, class TCleanupFunction = void(*)(TArgsInfo*) >
class args_info_manager
{
public:
	args_info_manager(TArgsInfo & args_info, TCleanupFunction cf)
	{
		this->args_info_pointer = &args_info;
		this->cleanup_function = cf;
	}
	~args_info_manager()
	{
		this->cleanup_function(this->args_info_pointer);
	}
private:
	TArgsInfo * args_info_pointer;
	TCleanupFunction cleanup_function;
};
//--------------------------------------------------------------------
/** \brief Process gengetopt with config file option (shamelessly stolen from RTK)
*
* \author Simon Rit
*
* \ingroup Macro
*/
#define GGO(ggo_filename, args_info)                                                                       \
	args_info_##ggo_filename args_info;                                                                    \
	cmdline_parser_##ggo_filename##_params args_params;                                                    \
	cmdline_parser_##ggo_filename##_params_init(&args_params);                                             \
	args_params.print_errors = 1;                                                                          \
	args_params.check_required = 0;                                                                        \
	args_params.override = 1;                                                                              \
	args_params.initialize = 1;                                                                            \
if (0 != cmdline_parser_##ggo_filename##_ext(argc, argv, &args_info, &args_params))                        \
	{                                                                                                      \
	std::cerr << "Error in cmdline_parser_" #ggo_filename "_ext" << std::endl;                             \
	exit(1);                                                                                               \
	}                                                                                                      \
	std::string configFile;                                                                                \
if (args_info.config_given)                                                                                \
	configFile = args_info.config_arg;                                                                     \
	cmdline_parser_##ggo_filename##_free(&args_info);                                                      \
if (configFile != "")                                                                                      \
	{                                                                                                      \
if (0 != cmdline_parser_##ggo_filename##_config_file(configFile.c_str(), &args_info, &args_params))        \
	  {                                                                                                    \
	  std::cerr << "Error in cmdline_parser_" #ggo_filename "_config_file" << std::endl;                   \
	  exit(1);                                                                                             \
	  }                                                                                                    \
	  args_params.initialize = 0;                                                                          \
	}                                                                                                      \
	args_params.check_required = 1;                                                                        \
if (0 != cmdline_parser_##ggo_filename##_ext(argc, argv, &args_info, &args_params))                        \
	{                                                                                                      \
	std::cerr << "Error in cmdline_parser_" #ggo_filename "_ext" << std::endl;                             \
	exit(1);                                                                                               \
	}                                                                                                      \
	args_info_manager< args_info_##ggo_filename >                                                          \
	manager_object(args_info, cmdline_parser_##ggo_filename##_free);
//--------------------------------------------------------------------