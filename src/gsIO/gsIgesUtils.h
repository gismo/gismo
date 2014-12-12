

#pragma once

#include <gsIO/gsXmlUtils.h>

namespace gismo {

namespace internal {

  // class gsXmlNode;
  // class gsXmlAttribute;
  // class gsXmlTree;

class gsIges
{
  typedef std::istream istream;
  typedef std::string string;
  struct IgesHeader;

private:
  typedef const char* ccp;
  gsIges(){ };

public:
  gsIges(istream& is, gsXmlTree * data);

  void readIges(istream& is);


	///Return directory entry point of DEid.
	MGIgesDirectoryEntry* directoryEntry(int DEid){
		return m_DirectoryEntries[DEid];
	};
	const MGIgesDirectoryEntry* directoryEntry(int DEid)const{
		return m_DirectoryEntries[DEid];
	};

	void clearStartSection(){m_StartSection=std::string();};
	void clearGSection(){m_GSection=MGIgesGSec();};
	void clearDirectoryEntries(){m_DirectoryEntries.clear();};
	void clear(){clearStartSection(); clearDirectoryEntries(); clearDirectoryEntries();};
	void set_initial_StartSection();

	const MGIgesGSec& GSection()const{return m_GSection;};
	MGIgesGSec& GSection(){return m_GSection;};

	///get the output line number of Start Section.
	int get_line_number_of_SS()const{return m_StartSection.size()/72+1;};

	///get the output line number of Global Sections.
	int get_line_number_of_GS()const{return m_nlineGSec;};

	///get the output line number of Directory Entries.
	int get_line_number_of_DE()const{return (m_DirectoryEntries.size()-1)*2;};

private: // Structs

  enum IgesSection { S, G, D, P, T, E };

/// Storage of all data contained in the IGES header 
  struct IgesHeader
  {
    IgesHeader() { };
    IgesHeader(string g);

    char m_delimeter_param;				// 1
    char m_delimeter_record;				// 2
    std::string m_productID_sender;			// 3
    std::string m_file_name;				// 4
    std::string m_native_systemID;			// 5
    std::string m_preprocessor_version;			// 6
    int m_number_of_bits_of_integer;			// 7
    int m_magnitude_single_precision;			// 8
    int m_significance_single_precision;		// 9
    int m_magnitude_double_precision;			//10
    int m_significance_double_precision;		//11
    std::string m_productID_receiver;			//12
    double m_model_space_scale;				//13
    int m_unit_flag;					//14
    std::string m_unit_name;				//15
    int m_max_number_of_line_weight_gradations;	        //16
    double m_width_of_max_line_weight;			//17
    std::string m_DateTime_File_generation;		//18
    double m_min_resolution;				//19
    double m_max_coordinate_value;			//20
    std::string m_author_name;				//21
    std::string m_author_organazation;			//22
    int m_version_flag;					//23
    int m_drafting_standard_flag;			//24
    std::string m_DateTime_Model_generation;	        //25
    std::string m_application_protocolID;		//26
  };

  /// Storage of all data contained in an IGES directory entity
  struct IgesDirectory
  {
    IgesDirectory(const char* start);

    int entity_type_number;
    int param_data_start;
    int structure;          // Typically 0
    int line_font_pattern;  // Typically 1 (solid)
    int level;              // Typically 0
    int view;               // Typically 0
    int trans_matrix;       // Typically 0
    int label_display;      // Typically 0
    std::string status;          // 8 digits, 
    double line_weight;        // Typically 0
    int color;              // Typically 0 (no color specified)
    int line_count;
    int form;               // 0 for free-form nurbs
    std::string entity_label;
    int entity_number;      // Typically 0
  };
  
  /// Supported entities (all others are ignored)
  bool supportedEntity(int ent)
    {
      return ( ent == 126 || ent == 128 || ent == 116 );
      // 100 Circular arc.
      // 102 Composite curve.
      // 104 Conic arc.
      // 106 Linear path (form 12).
      // 108 Plane.
      // 110 Line.
      // 116 Point.
      // 118 Ruled surface.
      // 120 Surface of revolution.
      // 122 Tabulated cylinder.
      // 123 Direction.
      // 124 Transformation matrix.
      // 126 Rational B-spline curve.
      // 128 Rational B-spline surface.
      // 141 Boundary.
      // 142 Curve on a parametric surface.
      // 143 Bounded surface.
      // 144 Trimmed (parametric) surface.
      // 186 Manifold solid b-rep object.
      // 190 Plane surface.
      // 192 Right circular cylindrical surface.
      // 314 Color definition.
      // 402 Group without back pointers assoc. (form 7).
      // 502 Vertex List.
      // 504 Edge List.
      // 508 Loop.
      // 510 Face.
      // 514 Shell.      
    }

private:

  void readIgesCurve(const char* start); 
  void readIgesSurface(const char* start);
  void readIgesPointCloud(const char* start);
    
  static int readIGESint(ccp& start, char pd, char rd); 
  static double readIgesDouble(ccp& start, char pd, char rd); 
  static string readIgesString(ccp& start, char delim, char delim2 = ';' );
  static bool readSingleIgesLine(istream& is, char line_terminated[81],
				 IgesSection& sect);

  static void skipDelimiter(ccp& whereami, char pd);
  static bool checkDelimiter(ccp& whereami, char wanted_delimiter,
		      char alternative_delimiter);
  static void skipOptionalTrailingArguments(ccp& whereami,
				     char pd, char rd,
				     int max_to_skip = 2); 
  
// ===========================================================

///Read in DE pointer into DEpointer.
///Line number in the istrm is converted to DE pointer.
///Function's return value is
///  true: when value specified.
///  false:when value not specified, DEpointer be 0.
bool get_DEpointer(
	char pDelimeter,	///<parameter delimeter
	std::istringstream& istrm,	///<Input string stream that contains integer data.
		///<The stream pointer will be advanced to the start position of the next item.
	int& DEpointer	///<output integer data that is converted from the istrm data.
);

///Read in Hollerith_string into strngData.
///Function's return value is
///  true: when value specified, strngData.size() be >0.
///  false:when value not specified, strngData.size() be 0.
bool get_Hollerith_string(
	char pDelimeter,	///<parameter delimeter
	std::istringstream& istrm,	///<Input string stream that contains Hollerith data.
		///<The stream pointer will be advanced to the start position of the next item.
	std::string& strngData///<output string data that is converted from the istrm's Hollerith data.
);

///Read in integer_string into intData.
///Function's return value is
///  true: when value specified.
///  false:when value not specified, intData be 0.
bool get_integer(
	char pDelimeter,	///<parameter delimeter
	std::istringstream& istrm,	///<Input string stream that contains integer data.
		///<The stream pointer will be advanced to the start position of the next item.
	int& intData	///<output integer data that is converted from the istrm data.
);
bool get_integer(
	char pDelimeter,	///<parameter delimeter
	std::istringstream& istrm,	///<Input string stream that contains integer data.
		///<The stream pointer will be advanced to the start position of the next item.
	short& shortData	///<output integer data that is converted from the istrm data.
);

///Read in real_string into realData
///Function's return value is
///  true: when value specified.
///  false:when value not specified, realData be 0.
bool get_real(
	char pDelimeter,	///<parameter delimeter
	std::istringstream& istrm,	///<Input string stream that contains real data.
		///<The stream pointer will be advanced to the start position of the next item.
	double& realData	///<converted real data from istrm will be output.
);
bool get_real(
	char pDelimeter,	///<parameter delimeter
	std::istringstream& istrm,	///<Input string stream that contains real data.
		///<The stream pointer will be advanced to the start position of the next item.
	float& floatData	///<converted real data from istrm will be output.
);

///convert the line id into int(sequence), inputting one line.
void get_ID_sequence(
	const std::string& line,///<Input whole line data(1-80)
	char& sectionID_letter,	///<section identification letter of the line will be output.
	int& sequence			///<ascending sequence number of the line.
);


private: // Data

  gsXmlTree * m_xmldata;
  std::ifstream m_ifstream;///<Input stream.

  std::string m_StartSection;///<Start section string data.
  IgesGlobal m_GSection;     ///<Global section data.
  int m_nlineGSec;	     ///<Number of the lines of Global section.
  std::vector<IgesDirEntry> m_DirectoryEntries;///<Directry entry data vector.
  ///<One pair of directory entry lines are stored in m_DirectryEntry[i].
				///<m_DirectoryEntries[0] is a dummy entry and has no meaning.

};


}// namespace internal

}// namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsIgesUtils.hpp)
#endif
