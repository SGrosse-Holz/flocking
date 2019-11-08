/* This is a small tool for easy command line option handling.
 * For a quick start, copy the following code into your
 * 	int main (int argc, char *argv[]) :

op::OptionHandler OH (R"HELP(
This is the help text (as C++11 raw string literal)
)HELP");
int istore = 0;
OH.addOption((new op::SingleValueOption<int>("i", istore))->description("Integer option"));
op::pRes res = OH.procOptions(argc, argv);
if (res != op::ok) return res;

 * This registers an option -i, which takes an integer argument and stores it in the variable istore.
 * Furthermore, the options -h[elp] and -? are always available; they display the help text.
 */

#include <map>
#include <string>
#include <iostream>
#include <sstream>

namespace op {

typedef enum {ok = 1, ret = 0, ret_err = -1} pRes;

class Option
{
	public:
				Option(std::string key, int count = 0)
					: argc_(count), desc_(""), name_(key) {}
				virtual ~Option() {}

		std::string	name() {return name_;}
		int		argc() {return argc_;}
		std::string	helpline()
				{
					if (desc_.length() > 0)
					{
						return std::string("  -").append(name_)
							.append("\t").append(desc_).append("\n\n");
					}
					return std::string("");
				}

		virtual
		Option*		description(std::string desc)
				{
					desc_ = desc;
					return this;
				}
		virtual
		pRes		proc(char *argv[]) = 0;

	protected:
		int		argc_;

	private:
		std::string	desc_;
		std::string	name_;
}; // class Option
		
class OptionHandler
{
	public:
			OptionHandler	(std::string helptext) : desc_(helptext) {}
		       ~OptionHandler	() {for (auto& x: opMap_) delete x.second;}
		void	addOption	(Option *opt) {opMap_[opt->name()] = opt;}
		pRes	procOptions	(int argc, char *argv[]);
		std::string	getHelp	()
				{
					std::string ret = "\n";
					ret += desc_;
					ret += "\n\nAvailable options:\n\n";
					ret += "  -h[elp]\n  -?\tDisplays this help\n\n";
					for (auto& x : opMap_)
						ret += x.second->helpline();
					return ret;
				}
	
	private:
		std::map<std::string, Option*>	opMap_;
		const std::string		desc_;
}; // class OptionHandler

template <class value_type>
class SingleValueOption : public Option
{
	public:
				SingleValueOption(std::string key, value_type& storage)
					: Option(key, 1), store_(storage) {}

		pRes		proc(char *argv[]);

	private:
		value_type&	store_;
}; // class SingleValueOption

template <class value_type>
pRes SingleValueOption<value_type>::proc (char *argv[])
{
	std::istringstream oss (argv[1]);
	oss >> store_;
	return ok;
}

pRes OptionHandler::procOptions (int argc, char *argv[])
{
	bool error = false;
	int errorAt = 0;
	for (int i = 1; i < argc; i++)
	{
		if (argv[i][0] != '-')
		{
			error = true;
			errorAt = i;
			break;
		}
		
		char *curName = argv[i]+1;
		if (	   std::string(curName) == "h"
			|| std::string(curName) == "help"
			|| std::string(curName) == "?"		)
		{
			std::cout << getHelp();
			return ret;
		}

		try
		{
			Option *opt = opMap_.at(argv[i]+1);
			if (i + opt->argc() >= argc)
			{
				std::cout << "Too few arguments specified for option " << argv[i]
					  << "; expected " << opt->argc() << std::endl;
				return ret_err;
			}
			pRes res = opt->proc(argv+i);
			if (res != ok)
				return res;

			i += opt->argc();
		}
		catch (const std::out_of_range &oor)
		{
			error = true;
			errorAt = i;
			break;
		}
	}

	if (error)
	{
		std::cout << "Invalid option: " << argv[errorAt] << std::endl;
		return ret_err;
	}
	return ok;
} // OptionHandler::procOptions

} // namespace op
