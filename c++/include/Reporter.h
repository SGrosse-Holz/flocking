/* This file is distributed under MIT license as part of the project flocksims.
 * See LICENSE file for details.
 * (c) Simon Grosse-Holz, 2019
 */

# ifndef REPORTER_H
# define REPORTER_H

# include "State.h"
# include "H5Cpp.h"
# include <string>
# include <list>
# include <iostream>

namespace flocksims {

namespace reportModes {
	enum reportMode {all=0, noParticles=1, noAnalysis=1<<1}; // The flags can be combined with bitwise or |
}
using reportModes::reportMode;

class BaseReporter
{
	public:
		BaseReporter(std::string filename);
		virtual ~BaseReporter() {};

		// Reporting functions
		template<typename T> // report name value pairs as attributes to root group
		void report(const char* name, const T& value, const H5::H5Object *parent=0);
		virtual void report(const State&, reportMode toReport=reportModes::all) = 0; // as individual data set

		// Write actual file
		virtual void dump() {outfile.close();}

	protected:
		H5::H5File outfile;
		std::string filename;
		bool do_report;

		// Overloaded function to resolve the template above
		H5::DataType DTfromValue(const int&) {return H5::DataType(H5::PredType::NATIVE_INT);}
		H5::DataType DTfromValue(const double&) {return H5::DataType(H5::PredType::NATIVE_DOUBLE);}
		H5::DataType DTfromValue(const std::string& val) {return H5::StrType(H5::PredType::C_S1, val.length());}
};

class HDF5Reporter : public BaseReporter
{
	public:
		HDF5Reporter (std::string filename) : BaseReporter(filename) {}

		// This inheritance has to be very explicit, otherwise
		// overloading and template will not play together
		template<typename T>
		void report(const char* name, const T& value, const H5::H5Object *parent=0)
		{ BaseReporter::report(name, value, parent); }

		void report(const State&, reportMode toReport=reportModes::all);

		void dump();

	private:
		std::list<State> datalist;
};

// template definition has to be in header, such that it is instantiated for
// use in other translation units. Example was void report<int>(...) which
// doesn't happen in HDF5Reporter.cpp but in main.cpp.
template<typename T>
void BaseReporter::report(const char* name, const T& val, const H5::H5Object *parent)
{
	outfile.openFile(filename, H5F_ACC_RDWR);

	if (!do_report) return;
	if (!parent) parent = &outfile;

	H5::DataType datatype = DTfromValue(val);
	H5::DataSpace dataspace(H5S_SCALAR);
	H5::Attribute attribute = parent->createAttribute(name, datatype, dataspace);
	if (std::is_same<T, std::string>::value)
		attribute.write(datatype, ((std::string*)&val)->c_str());
	else
		attribute.write(datatype, &val);

	outfile.close();
}

}

# endif
