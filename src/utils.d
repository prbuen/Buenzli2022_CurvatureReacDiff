/*
  Utils: dumping place for some ugly things written through the years
  
  Author: Pascal R Buenzli, 2022
*/

module utils;

import std.stdio;
import std.conv;
import std.exception;
import std.traits; // : has_member
// static import std.traits;
import std.regex;
import std.path;
import std.algorithm;
import std.math;
import std.process;
import std.string;
import std.file; // : exists, isFile;
import std.typecons;
import std.range;
import std.format;


/*
To select a version and export this choice, do
---
version=Name
mixin(DefSwitch!("Name"))
---
The mixin allows the enum bool `Name` to be used as a (`static if`) switch in other modules that import this one. This mixin is a shorthand for writing
---
version(Name)
	enum bool Name = true;
else
	enum bool Name = false;
---
 */
template DefSwitch(string Name)
{
	const char[] DefSwitch = "version(" ~ Name ~ ") enum bool " ~ Name ~ " = true; else enum bool " ~ Name ~ " = false;";
}

/*
Emulation of a string output stream

Example:
-----
sstring s = "Hello";
s ~= ", ";
int n = 10;
s.writeln("world, ", n, " times");
s.writefln(" and %d times", n);
-----
 */
struct sstring {
	string _s;
	alias _s this; // treat sstring just like _s, a string
	
	this(string literal) { _s = literal; }

	void write(Args...)(Args args) { foreach(a; args)  {_s ~= to!string(a);} }
	void writeln(Args...)(Args args) { this.write(args, '\n'); }
	void writef(Char, Args...)(in Char[] fmt, Args args) { _s ~= std.string.format(fmt, args); }
	void writefln(Char, Args...)(in Char[] fmt, Args args) { this.writef(fmt ~ '\n', args); }	
}


/*_________________________________________[ Data ]_________________________________________ */
/**
Convenience structure to simplify opening files for writing the state. Contains information about the name of the directory containing the datafiles (required minimum), and (optionally) a pattern for datafile names, such as <basename>, digital format, and <ext>.

For a datafile in datadir called 'mydata_000.dat', we disentangle name_pattern = 'mydata', digit_pattern='3', ext_pattern='.dat', number_format = "_%03d" (integers left-padded with 0 for a total width of 3). If the datafile is simply mydata_0.dat, then digit_pattern='0' and number_format="%d". This is pure pattern parsing, no check is done for the existence of files/directory.
 */
struct Data {
public:
	string datadir; // directory containing datafiles
	string name_pattern; // basename, i.e., filename witout _<frame_id>.<ext>
	size_t digit_pattern; // number of digits if frame_id of the form _0001 (=4), 0 if of the form _0, ..., _9, _10, etc.
	string ext_pattern; // file extension
	string number_format; // so we can use sstring s; s.writef(number_format, frame_id) to encode the frame_id in the right format (prepending with 0 if this is how digit_pattern is initialised)
	size_t init_frame_id; // initial frame id. Usually 0, but can be larger (positive)

	// initialise datadir, name_pattern, digit_pattern, ext_pattern, and number_format by analysing the filename pattern: <name>_0.<ext> or <name>_00.<ext> etc.
	this(in string firstdatafile)
	{
		string filename = baseName(firstdatafile);
		datadir = dirName(firstdatafile);
		ext_pattern = extension(filename);

		auto m = matchFirst(filename, regex(`(?P<name>[\.\w\s]*)_(?P<init_number>[0-9]*)\.(?P<ext>[\w\s\.]*)`));

		enforce(!m.empty, std.conv.text("*** Data class initialisation. Invalid frame data filename: ", filename, ". The filename must have format '<name>_<number>.<ext>'"));

		name_pattern = m["name"];
		string numberstring = m["init_number"];
		init_frame_id = to!size_t(numberstring);
		digit_pattern = 0;
		if (numberstring.length > 1 && numberstring[0]=='0')
			digit_pattern = numberstring.length;

		sstring number_format_tmp = "";
		if (digit_pattern > 0)
		{
			number_format_tmp = `_%0`;
			number_format_tmp.writef("%dd", digit_pattern);
		}
		else
			number_format_tmp = `_%d`;
		number_format = number_format_tmp;

		ext_pattern = "." ~ m["ext"];
		
		// //###############
		// stderr.writeln("Data directory: ", datadir);
		// stderr.writeln("Name pattern: ", name_pattern);
		// stderr.writeln("Initial frame no.: ", init_frame_id);
		// stderr.writeln("No. of digits (0 if variable): ", digit_pattern);
		// stderr.writeln("Number format: ", number_format);
		// stderr.writeln("ext: ", ext_pattern);
		// enforce(0, "FINISHING HERE TO TEST"); //#######
		// //###############
	}

	this(in string datadir_, in string name_pattern_="", in string ext_pattern_="", in string numberstring="0")
	{
		datadir = datadir_;
		name_pattern = name_pattern_;
		ext_pattern = ext_pattern_;

		init_frame_id = to!size_t(numberstring);
		digit_pattern = 0;
		if (numberstring.length > 1 && numberstring[0]=='0')
			digit_pattern = numberstring.length;

		sstring number_format_tmp = "";
		if (digit_pattern > 0)
		{
			number_format_tmp = `_%0`;
			number_format_tmp.writef("%dd", digit_pattern);
		}
		else
			number_format_tmp = `_%d`;
		number_format = number_format_tmp;		
	}

	/* ===== data.filepath() =====
    /*
	use as
	    filename = data.filepath('data', 10, '.dat') => returns the path '<datadir>/data_010.dat', appropriately normalised (including \ instead of / on windows)
	*/
	string filepath(in string name, in size_t frame, in string ext_="")
	{
		string ext = ext_;
		if (ext=="")
			ext = ext_pattern;
		return buildNormalizedPath(datadir, name ~ std.string.format(number_format, frame) ~ ext);
	}

    /*
	use as
	    filename = data.filepath('metadata.dat') => returns the path '<datadir>/metadata.dat', appropriately normalised (including \ instead of / on windows)
	*/
	string filepath(in string filename)
	{
		return buildNormalizedPath(datadir, filename);
	}
	
	/* ===== data.open() =====
	/*
	use as
		auto data_p = data.open('data', 10, '.dat'); // opens '<datadir>/data_010.dat' for writing in binary mode (to prevent \n -> \r\n conversion in windows)
		data_p.write(...);
	 */
	File open(in string name, in size_t frame, in string ext_="", string rwmode="wb")
	{
		string ext = ext_;
		if (ext=="")
			ext = ext_pattern;
		return File(buildNormalizedPath(datadir, name ~ std.string.format(number_format, frame) ~ ext), rwmode);
	}

	/* ===== data.open() =====
	/*
	use as
		auto metadata_p = data.open('metadata.dat'); // opens '<datadir>/metadata.dat' for writing in binary mode (to prevent \n -> \r\n conversion in windows)
		metadata_p.write(...);
	 */
	File open(in string filename, string rwmode="wb")
	{
		return File(buildNormalizedPath(datadir, filename), rwmode);
	}


	/*===== erase_datadir() =====*/
	/*
	   Erases the datafiles directory, after a (nonexhaustive) number of security checks on data.datadir. This can be called by the user in the Model class constructor.

	   The datafiles directory `datadir` can't be '.', can't be '/' (they would be erased!), and it can't be a symlink (because it would be erased and replaced with a real directory).

	   If `forbid_parent` is `true`, `datdir` cannot start with '..' (even if not actually referring to a parent, but e.g. branching off to an independent path). If `forbid_absolute` is `true`, `datadir` cannot be an absolute path.

	   It is recommended that `forbid_parent` and `forbid_absolute` are maintained at their default value of `true`, to enforce that the datafiles directory is a real subdirectory of the executable's directory.
	 */
	void erase_datadir(bool forbid_parent = true, bool forbid_absolute = true)
	{
		// Safety:
		enforce(buildNormalizedPath(datadir) != buildNormalizedPath("."), "Cannot choose '.' for the datafiles directory, see erase_datadir().");
		version(Posix)
			enforce(buildNormalizedPath(datadir) != "/", "Cannot choose root '/' for the datafiles directory, see erase_datadir().");
		version(Windows)
			if(buildNormalizedPath(datadir).length == 1)
				enforce(! isRooted(buildNormalizedPath(datadir)), "Cannot choose root '\' for the datafiles directory, see erase_datadir().");
			else if(buildNormalizedPath(datadir).length <= 3)
				enforce(! isRooted(buildNormalizedPath(datadir)), std.conv.text("Cannot choose root '", datadir, "' for the datafiles directory, see erase_datadir()."));
		if (exists(datadir))
			enforce(!isSymlink(buildNormalizedPath(datadir)), std.conv.text("Cannot use a symlink ('", datadir, "') for the datafiles directory (it would be erased and replaced with a real directory)."));

		if (forbid_parent && datadir.length > 1)
			enforce(buildNormalizedPath(datadir)[0 .. 2] != "..", std.conv.text("Cannot choose parent directory (or subdirectory thereof) for the datafiles directory ('", datadir, "'), see erase_datadir(). (Use symbolic link to executable in parent directory of datafiles directory instead.)"));

		if (forbid_absolute)
			enforce(! isAbsolute(datadir) && ! isRooted(datadir), std.conv.text("Cannot use an absolute path ('", datadir, "') for the datafiles directory, see erase_datadir()."));			
		
		// Erase all existing data files in datadir and (re)create 
		if(exists(datadir))
		{
			// stderr.writeln("Removing directory '", datadir, "'"); //#########
			// enforce(0); //##########
			rmdirRecurse(datadir);
			enforce(! exists(datadir), std.conv.text("Could not remove data directory '", datadir, "'")); // sometimes rmdir fails if too many data files?
		}
		mkdirRecurse(datadir);
	}	
}

enum string[1] All = [""];
enum string[] Functions = [" All functions "]; // any name we can test against and that can't be an identifier.
enum string[] None = [];

/*____________________________________[ ModelInputsWrap ]______________________________________*/
/**
   Classes implementing the `ModelInputs` base class will define the parameters and input functions of the computer model. The usefulness of the wrapper ModelInputsWrap class is in adding `read_params()` to read parameter values from a config file, and params() to write current parameter values in a formatted way into a data file. This wrapper is necessary so that traits are still able to introspect the fields of the ModelInputs class at compile time. Inheriting from a base ModelInputs class that would define params() and read_params() doesn't allow introspection of the children, because inheritance is run-time, but traits are compile time. However, the 'alias this' enables to mimick this behaviour well.

   If `ModelInputsT` declares a string `metadata`, a metadata file `datadir/metadata.dat` will automatically be written with `p.metadata`'s content. Otherwise, the data file directory will contain no file `metadata.dat`. The user defines the content of `p.metadata`. Typically, it contains a list of parameter values used to run the simulation. The method `read_params()` can be used to read parameter values from text files to `p`, and the method `params()` can be used to write model parameter values to the string `p.metadata` during initialisation, prior to it being written to the metadata file.

   If `ModelInputsT` defines the string `model_version` and `Model`'s constructor is given a string argument (e.g. the executable name), it is checked whether `model_version` is contained within that string. If not, a warning is issued on `stderr`. This is useful to check potential mismatches between the file structure (which may contain a version number) and `model_version` defined in the source code (which may be outputted to metadata, copyright notices, etc).
*/


class ModelInputsWrap(ModelInputsT) {
  public:
	ModelInputsT p;
	alias p this;
  public:
	/*===== this() =====*/
	this()
	{
		static if ( is (ModelInputsT == struct) ) // http://dlang.org/expression.html
			p = ModelInputsT();
		else
			p = new ModelInputsT;
	}
	
	/*===== params() =====*/
	/** params(comment_prefix, commented_list, uncommented_list, title, count)
	   Returns a string listing the model specifications (name of derived Model class, and name of template class parameter ModelInputsT) and then listing all the parameter values of the model contained in the instance `p` (except the special variable `p.metadata`). This is a convenience function that can be used e.g. in `Model.ctor()` to print out parameter values associated with the model's data output. The special variable `p.metadata` is not included in the output of `params()`. This variable is usually set to the output of `params()` by the user, so excluding `p.metadata` from the output of `params()` prevents bloating the output upon successive calls to `params()`.  
	*/
	/*(This function is deferred to `Model` rather than provided in a `ModelInputs` interface because nonstatic variables need `ModelInputs` to be instantiated first.)*/
	/**
	   Params:
			comment_prefix = string with which to comment out writing a variable and its numeric value (scalars and arrays of scalars) or callables (functions) and other nonnumeric quantities (class instances)
			commented_list = list of items to prepend with `comment_prefix`
			uncommented_list = list of items $(B not) to prepend with `comment_prefix` (`uncommented_list` has precedence over `commented_list`)
	   		title = whether to 'pretty print' with titles/subtitles
			count = whether to print a summary of parameter count at the end

		Special_lists: For `commented_list` and `uncommented_list`, one can use 'None' to mean an empty list, 'All' to mean all the possible items, and 'Functions' to mean all the callable items. If `comment_prefix` is changed to `""` these lists have no effect, nothing can be commented out (not recommended). $(BR)
	   
		Examples:
		The most common use cases will be:
		----
		params(); // nothing commented out
		params("# "); // all commented out with "# "
		params("# ", ["description", "string1", "string2"]); // only `description`, `string1`, `string2` commented out with "# "
		params("# ", ["description", "string1"] ~ Functions); // `description`, `string`, and all functions commented out with "# "
		params("# ", All, ["dt", "dx"]) // all commented out with "# " except `dt` and `dx` (uncommented_list > commented_list)
		----

		Further examples, listing alternative syntax:
		-----
		// Nothing commented out:
		params()
		params("# ", None); // Title and count are preceded with "# "
		params("# ", None, All);
		params("# ", None, None);
		params("# ", All, All); // uncommented_list > commented_list for overlaps
		params("# ", ["param1", "param2"], All);

		// Some commented out with "# ": `param1` and `param2`, the others not commented out
		params("# ", ["param1", "param2"]);
		params("# ", ["param1", "param2"], None);
		params("# ", All, ["param3", "param4", "param5"]); // assuming there is `param1`, `param2`, `param3`, `param4`, `param5`
		params("# ", ["param1", "param2", "param3"], ["param3", "param4", "param5"]); // assuming there is `param1`, `param2`, `param3`, `param4`, `param5`

		// Functions commented out
		params("# ", Functions);
		params("# ", ["param1", "param2"] ~ Functions);

		// All commented out except `param4`, `param5`:
		params("# ", All, ["param4", "param5"]);

		// All commented out:
		params("# ");
		params("# ", All);
		params("# ", All, None);
		-----

		Note_for_Windows_users: Parameters are separated by unix-like end-of-line characters `'\n'` even under Windows. An automatic conversion of `'\n'` to `'\r\n'` may occur when writing the output of `params()` to a file/stream open in text mode. See `pb.windows.setFileModeBinary` to force writing in binary mode and prevent this conversion.
*/
	final string params(in string comment_prefix = "", in string[] commented_list = All, in string[] uncommented_list = None, in bool title = true, in bool count = true)
	{

		string prefix_title = comment_prefix; // use the same comment prefix for title and count		
		// stderr.writeln("commented_list = ", commented_list);
		// stderr.writeln("uncommented_list = ", uncommented_list);
		
		// helper functions:
		string pre_line(in string comment_prefix, in string name, in bool is_function = false)
		{
			if (commented_list == None)
				return "";
			
			if (uncommented_list == All || (uncommented_list != None && uncommented_list.canFind(name))) // name in uncommented_list (has precedence over name in commented_list, so to be tested first)
				return "";

			if (is_function && uncommented_list != None && uncommented_list.canFind(Functions[0]))
				return "";

			if (commented_list == All || commented_list.canFind(name) || (is_function && commented_list.canFind(Functions[0])))
				return comment_prefix;
			else
				return "";			
		}
						
		// Writing the name and value of all relevant fields of ModelInputsT:
		// see https://issues.dlang.org/show_bug.cgi?id=12791 for dealing with inherited members
		// see http://forum.dlang.org/thread/zubtrrbmphjcapgheljl@forum.dlang.org for dealing with private members
		string s="";

		// prepend `metadata` by the name of the `ModelInputs` class so this is recorded in the datafile automatically
		if(title) s ~= prefix_title ~ "===== Model =====" ~ "\n";
		s ~= comment_prefix ~ std.conv.text(typeid(ModelInputsT)) ~ "\n";
		s ~= comment_prefix ~ std.conv.text(typeid(this)) ~ "\n";
		if (title) s ~= comment_prefix ~ "\n";
		
		if (title) s ~= prefix_title ~ "===== Model inputs =====\n";
		
		size_t count_var=0;
		size_t count_fun = 0;
		
		foreach(i, m; __traits(allMembers, ModelInputsT))
		{
			// stderr.writeln("m=", m);//#####################
			
			static if (__traits(compiles, std.conv.text(__traits(getMember, p, m)))
					   && m != "__ctor"
					   && m != "toString"
					   && m != "toHash"
					   && m != "params"
					   && m != "read_params" // getMember of p, not of ModelInputsT to see nonstatic members!
					   && ! isCallable!(__traits(getMember, ModelInputsT, m)) // functions that return dynamic arrays compile with the above but are to be dealt with as callables below, so we prevent that case
				)
			{
				// stderr.writeln("typeof(m)=", typeof(__traits(getMember, p, m)).stringof);//#####################
				static if(typeof(__traits(getMember, p, m)).stringof == "string" || typeof(__traits(getMember, p, m)).stringof == "char" || typeof(__traits(getMember, p, m)).stringof == "immutable(string)" ||
				typeof(__traits(getMember, p, m)).stringof == "sstring")
				{
					static if(m != "metadata")
						s ~= pre_line(comment_prefix, m) ~ m ~ " = \"\"\"" ~ __traits(getMember, p, m) ~ "\"\"\"" /* ~ "(of type " ~ typeof(__traits(getMember, p, m)).stringof ~ ")" */ ~ "\n";
				}
				else
				{
					// pragma(msg, "Debug: m=", m.stringof, " is of type: ", typeof(__traits(getMember, p, m)).stringof); // ###################

					s ~= pre_line(comment_prefix, m) ~ m ~ " = " ~ std.conv.text(__traits(getMember, p, m)) /* ~ "(of type " ~ typeof(__traits(getMember, p, m)).stringof ~ ")" */  ~ "\n";
					static if(typeof(__traits(getMember, p, m)).stringof != "string")
					{
						// stderr.writeln(m, " = ", __traits(getMember, p, m), " == ", typeof(__traits(getMember, p, m)).nan, "? ", std.conv.text(__traits(getMember, p, m)) == "nan");
						// comparing __traits(getMember,p,m) == typeof(__traits(getMember,p,m)).nan doesn't work...
						enforce("nan" != std.conv.text(__traits(getMember, p, m)), std.conv.text("***** ", m, " = ", __traits(getMember, p, m), " is not a number *****")); // throw exception if a parameter is NaN. TODO: Doesn't work for complex numeric types which would need to be tested with float.nan + i*float.nan
					}
				}
				count_var++;
			}
			else static if (__traits(compiles, __traits(getMember, ModelInputsT, m)) && m != "Monitor" && m != "factory" && m != "opCmp" && m != "opEquals" && m != "toHash" && m != "toString" && m != "__ctor") // first part is skip errors when m is a private member
			{
				static if ( isCallable!(__traits(getMember, ModelInputsT, m)) )
				{
					// stderr.writeln("########### m=", m, " IS CALLABLE ##############");//###########					
					s ~= pre_line(comment_prefix, m, true) ~ m ~ "(";
					alias ArgTypes = ParameterTypeTuple!(__traits(getMember, ModelInputsT, m));
					alias ArgIds = ParameterIdentifierTuple!(__traits(getMember, ModelInputsT, m));
					alias RetType = Unqual!(ReturnType!(__traits(getMember, ModelInputsT, m)));
					static assert (ArgTypes.length == ArgIds.length);
					foreach(var; ArgIds)
						s ~= var ~ ",";
					static if (ArgIds.length > 0)
						s = s[0 .. $-1];
					s ~= ") = ";
					static if (is (RetType == immutable(char)[]))
						s ~= "string(";
					else
						s ~= std.conv.text(typeid(RetType)) ~ "(";
					foreach(k, type; ArgTypes)
						static if (is (type == immutable(char)[]))
							s ~= "string" ~ " " ~ ArgIds[k] ~ ", ";
						else
							s ~= std.conv.text(typeid(Unqual!(type))) ~ " " ~ ArgIds[k] ~ ", ";
					static if (ArgTypes.length > 0)
						s = s[0 .. $-2];
					s ~= ")\n";// ~ " = " ~ std.conv.text(typeid(typeof(__traits(getMember, ModelInputsT, m)))) ~ "\n";
					count_fun++;
				}
				else
				{
					// pragma(msg, m, ", isAbstractFunction? ", isAbstractFunction!(__traits(getMember, ModelInputsT, m)));//##########
					// pragma(msg, m, ", isDelegate? ", isDelegate!(__traits(getMember, ModelInputsT, m)));//##########

					// pragma(msg, m, " getOverloads: ", __traits(getOverloads, ModelInputsT, m));
					// pragma(msg, m, " isSomeFunction? ", isSomeFunction!(__traits(getMember, ModelInputsT, m)));//###########
					// pragma(msg, m, " isTypeTuple? ", isTypeTuple!(__traits(getMember, ModelInputsT, m)));//###########					
					// pragma(msg, m, " isFinalFunction? ", isFinalFunction!(__traits(getMember, ModelInputsT, m)));//###########
					// pragma(msg, m, " isFunction? ", isFunction!(__traits(getMember, ModelInputsT, m)));//###########
					// pragma(msg, m, " compiles? ", __traits(compiles, {__traits(getMember, ModelInputsT, m)("a","b");}));

					s ~= pre_line(comment_prefix, m, true) ~ m ~ " = " ~ " <no value> and not callable, or templated delegate." ~ "\n";
				}
			}
		}

		
		if (count) s ~= prefix_title ~ std.conv.text("_____\n", prefix_title, count_var, " parameter(s) (array elements not included)\n", prefix_title, count_fun, " function(s)\n");
		// return s[0 .. $-1]; // remove trailing \n
		return s;
	}




	/*===== read_params() =====*/
	/**
	   Convenience function to import parameter values from a plain text file into fields of the variable instance `p` of type `ModelInputsT`. Calling `read_params!preamble(config_filename, names_to_ignore)` is equivalent to `p.read_config_mixin!preamble(config_filename, names_to_ignore ~ "metadata")` from `pb.config` if `p` contains the `mixin read_config_mixin;` template declaration (recommended due to scoping of imports). Otherwise, equivalent to `p.read_config!preamble(config_filename, names_to_ignore ~ "metadata")`. In the first case, most likely the preamble can be left blank, whereas it might need import declarations in the second case, see pb.config.read_config.

	   Note: `p.metadata` is always ignored - use `p.read_config[_mixin](...)` instead if you need to read `p.metadata` from a file).

	   It is possible to call the function either specifying a filename, or specifying a file pointer (this is useful to read from stdin):
	   ---
	   read_params("test.conf");
	   read_params(File("test.conf"));
	   read_params(stdin);
	   ---
	*/
	final void read_params(string preamble="")(in string config_filename, in string[] names_to_ignore = [])
	{
		// pragma(msg, "Is pa.read_config_mixin defined (filename)? ", __traits(compiles, {p.read_config_text_mixin("");}));//###############
		static if(__traits(compiles, {p.read_config_mixin("");}))
			p.read_config_mixin!preamble(config_filename, names_to_ignore ~ "metadata");
		else
			p.read_config!preamble(config_filename, names_to_ignore ~ "metadata");
	}
	final void read_params(string preamble="")(File config_file, in string[] names_to_ignore = [])
	{
		// pragma(msg, "Is p.read_config_mixin defined (file)? ", __traits(compiles, {p.read_config_mixin(config_file, names_to_ignore);}));//###############
		static if(__traits(compiles, {p.read_config_mixin(config_file, names_to_ignore);}))
			p.read_config_mixin!preamble(config_file, names_to_ignore ~ "metadata");
		else
			p.read_config!preamble(config_file, names_to_ignore ~ "metadata");
	}	
}


/**
  Vector:
  Simple one-dimensional array (implemented as dynamic 1D array to allow for large arrays contiguous in memory). No particular facility to deal with boundaries.
  Example:
  ---
  auto vec = new Vector!(100, double);
  vec[0] = ...
  ---

 */
class Vector(ValueType) {
	ValueType[] _M; // one dimensional array stored as _dynamic_ (otherwise only small arrays possible). Contiguous in memory
	string col_sep="\t"; // see toString

	alias _M this;
	this(in size_t N) { _M = new ValueType[N]; }

	/**
	   Outputting function allowing for writing as formatted string or raw output. See $(LINK http://wiki.dlang.org/Defining_custom_print_format_specifiers)$(BR)
	   Example:
	   ---
	   Vector!(Nx, double) m;
	   // initialise m...
	   writeln(m); // same as writefln("%s", m)
	   writefln("%s", m); // ascii output, each element printed
	   writefln("%r", m); // raw output, each element printed
	   writefln("%4.1s", m); // prints only elements at indices i*4+1, i=0,1,2,...
	   writefln("%4.1r", m); // prints only elements at indices i*4+1, i=0,1,2,...
		---
	*/
	void toString(scope void delegate(const(char)[]) sink, FormatSpec!char fmt) const
	{
		if (fmt.precision == fmt.UNSPECIFIED)
			fmt.precision = 0;
		if (fmt.width == 0)
			fmt.width = 1;
		
		if (fmt.spec == 'r') // raw format
			foreach(i; iota(fmt.precision, _M.length, fmt.width)) // i=fmt.precision, fmt.precision + fmt.width, fmt.precision + 2*fmt.width, ... < _M.length
				formattedWrite(sink, "%r", _M[i]);
		else if (fmt.spec == 's') // string (ascii) format
			foreach(i; iota(fmt.precision, _M.length, fmt.width))
				sink(std.conv.text(_M[i], col_sep));
		// {
		// 	// A col_sep is appended at the end. If this is a problem, do like this:
		// 	auto buffer = appender!string(); // preferrable to string buffer="";
		// 	foreach(i; iota(fmt.precision, _M.length, fmt.width)) // i=fmt.precision, fmt.precision + fmt.width, fmt.precision + 2*fmt.width, ... < _M.length
		// 		buffer ~= std.conv.text(_M[i], col_sep);
		// 	sink(buffer.data[0 .. $-1]);
		// }
		else
			throw new Exception("Unknown format specifier: '%" ~ fmt.spec ~ "'");
	}
	
}



/*___________________________[ adjust_to_dx_coarse_given_dx ]_____________________________*/
/**
   Adjust domain discretisation to `dx_coarse` given `dx` (both `domain_length` and `dx_coarse` will be adjusted to multiples of `dx`) and calculate related dependent parameters. This scheme is appropriate for temporal axis discretisation for which the end time point value is less important than the size of the time step.

   Params:
   domain_length = wished length of the domain; will be rounded up to a multiple of `dx_coarse` (in and out)
   dx = precision step size (in)
   Nx = total number of points that fit within `domain_length` with regular spacing `dx`. The first and last points coincide with the domain boundaries (out)
   dx_coarse = coarse step size between data output, may be adjusted to fit a multiple of `dx` (in and out)
   Nx_coarse = total number of points that fit within `domain_length` with regular spacing `dx_coarse`. The first and last coarse points coincide with the domain boundaries (out)
   output_every_nth_step = determines the frequency of data output across the domain (integer `dx_coarse/dx`) (out)
*/	
void adjust_to_dx_coarse_given_dx(Real) (ref Real domain_length, in Real dx, ref size_t Nx, ref Real dx_coarse, ref size_t Nx_coarse, ref size_t output_every_nth_step)
{
		output_every_nth_step = to!size_t(round(dx_coarse/dx));
		
		// re-adjust dx_coarse to align it with a multiple of dx:
		if (abs(dx_coarse/(output_every_nth_step*dx) - 1) > 0.01)
			stderr.writeln("***** WARNING *****\ndx_coarse = ", dx_coarse, " was changed to dx_coarse = ", dx*output_every_nth_step, " to align it with a multiple of dx = ", dx, ".\n This is more than 1% change on the initial value\n*******************\n");
		dx_coarse = output_every_nth_step * dx;

		// re-adjust domain_length to align it with a multiple of dx_coarse:
		Nx_coarse = 1 + to!size_t(ceil(domain_length/dx_coarse)); // // add 1: # points = # intervals + 1
		
		if (abs((Nx_coarse-1)*dx_coarse/domain_length - 1) > 0.01)
			stderr.writeln("***** WARNING *****\ndomain_length = ", domain_length, " was changed to domain_length = ", (Nx_coarse-1)*dx_coarse, " to align it with a multiple of dx*output_every_nth_step = ", dx_coarse, " (rounding up).\nThis is more than 1% change on the initial value\n*******************\n");
		domain_length = (Nx_coarse-1)*dx_coarse;
		Nx = 1 + to!size_t(ceil(domain_length/dx)); // add 1: # points = # intervals + 1
}


/*___________________________[ adjust_to_dx_coarse_given_length ]_____________________________*/
/**
   Adjust domain discretisation to `dx_coarse` given the domain length. This scheme is appropriate for spatial discretisation in which a given domain length is usually considered (as opposed to the time axis). First `dx_coarse` is adjusted so that it fits an integer number of times within domain_length (see `Nx_coarse`). Then `dx` is adjusted so that it fits an integer number of times within `dx_coarse` (see `Nx`).

   Params:
   domain_length = fixed length of the domain (in) = `x_max` - `x_xmin`
   dx = precision step size; its initial value will be adjusted (up or down) to fit an integer number of times within `dx_coarse`, see `Nx` (in and out)
   Nx = total number of points that fit within `domain_length` with regular spacing `dx`. The first and last points coincide with the domain boundaries (out). The number of intervals is Nx-1.
   dx_coarse = coarse step size between data output; its initial value may be adjusted (up or down) to fit an integer number of times within domain_length, see `Nx_coarse` (in and out)
   Nx_coarse = total number of points that fit within `domain_length` with regular spacing `dx_coarse`. The first and last points coincide with the domain boundaries (out). The number of intervals is Nx_coarse-1.
   output_every_nth_step = determines the frequency of data output across the domain (integer `dx_coarse/dx`) (out)
*/	
void adjust_to_dx_coarse_given_length(Real) (in Real domain_length, ref Real dx, ref size_t Nx, ref Real dx_coarse, ref size_t Nx_coarse, ref size_t output_every_nth_step)
{
	Nx_coarse = 1 + to!size_t(round(domain_length/dx_coarse)); // add 1: # points = # intervals + 1
	if (abs((Nx_coarse-1)*dx_coarse/domain_length - 1) > 0.01)
		stderr.writeln("***** WARNING *****\ndx_coarse = ", dx_coarse, " was changed to  = ", domain_length/(Nx_coarse-1), " so it fits an integer number of times within domain_length = ", domain_length, ".\nThis is more than 1% change on the initial value\n*******************\n");
	dx_coarse = domain_length/(Nx_coarse - 1); // adjust
	
	Nx = 1 + (Nx_coarse-1)*to!size_t(round(dx_coarse/dx)); // add 1: # points = # intervals + 1
	if (abs((Nx-1)*dx/domain_length - 1) > 0.01)
		stderr.writeln("***** WARNING *****\ndx = ", dx, " was changed to  = ", domain_length/(Nx-1), " so it fits an integer number of times within dx_coarse = ", dx_coarse, ".\nThis is more than 1% change on the initial value\n*******************\n");
	dx = domain_length/(Nx-1);

	output_every_nth_step = to!size_t(round(dx_coarse/dx));
}


/*
  Returns grad(phi) and Hess(phi) using 2nd order centred FD assuming periodic boundary conditions. The results returned in `Dx_phi` and `Dy_phi` etc. are NOT divided by the lattice space steps (dx, dy).
  Params: `phi` is the level set function, a `scid` real matrix with `phi.rows` rows and `phi.cols` columns $(BR)
  i,j the location in space where the approximation is sought
  Returns:
  `Dx_phi` = ( phi(i+1,j) - phi(i-1,j) )/2    (with PBC. Result is not divided by dx)
  `Dy_phi` = ( phi(i,j+1) - phi(i, j-1) )/2    (with PBC. Result is not divided by dx)
  `Dxx_phi` = phi(i+1, j) - 2phi(i,j) + phi(i-1,j)   (with PBC. Result is not divided by dx^2)
  `Dxy_phi` = Dx_centred (Dy_centred phi) = (1/4)(phi(i+1, j+1) - phi(i-1,j+1) - phi(i+1,j-1) + phi(i-1,j-1)) (with PBC. Result is not divided by dx dy)
  `Dyy_phi` = phi(i, j+1) - 2phi(i,j) + phi(i,j+1)   (with PBC. Result is not divided by dy^2)
*/
void grad_hess_phi_centred_2(MaterialProp,Real)(ref MaterialProp phi, in size_t i, in size_t j, ref Real Dx_phi, ref Real Dy_phi, ref Real Dxx_phi, ref Real Dxy_phi, ref Real Dyy_phi)
{
	auto Nx = phi.rows;
	auto Ny = phi.cols;

	// PBC:
	auto im1 = (i==0? Nx-1 : i-1);
	auto ip1 = (i==Nx-1? 0 : i+1);

	auto jm1 = (j==0? Nx-1 : j-1);
	auto jp1 = (j==Nx-1? 0 : j+1);

	// Calculate grad phi
	Dx_phi = 0.5*( phi[ip1,j] - phi[im1,j] );
	Dy_phi = 0.5*( phi[i, jp1] - phi[i, jm1] );

	// Calculate phi_xx, phi_xy, phi_yx, phi_yy
	Dxx_phi = phi[ip1,j] - 2*phi[i,j] + phi[im1,j];
	Dyy_phi = phi[i,jp1] - 2*phi[i,j] + phi[i,jm1];
	Dxy_phi = 0.25*( phi[ip1,jp1] - phi[im1,jp1] - phi[ip1,jm1] + phi[im1,jm1] );
}

/*
  Returns the signed curvature div(grad phi/| grad phi |) using centred 2nd order finite difference for all derivatives involved. See Sethian (1999), Eq (6.35) p69 for the formula. The result returned is NOT divided by the lattice step (dx or dy). This assumes dx=dy.
  Arguments: `Dx_phi` and `Dy_phi` returned are 2nd order centred finite difference (with PBC)
*/
Real curvature_centred_2(MaterialProp,Real)(ref MaterialProp phi, in size_t i, in size_t j, ref Real Dx_phi, ref Real Dy_phi)
{
	Real Dxx_phi, Dxy_phi, Dyy_phi;
	grad_hess_phi_centred_2(phi, i, j, Dx_phi, Dy_phi, Dxx_phi, Dxy_phi, Dyy_phi);
		
	return (Dxx_phi*Dy_phi^^2 - 2*Dy_phi*Dx_phi*Dxy_phi + Dyy_phi*Dx_phi^^2) / sqrt(Dx_phi^^2+Dy_phi^^2)^^3;
}
