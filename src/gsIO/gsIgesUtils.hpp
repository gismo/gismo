
#pragma once



namespace gismo {

namespace internal {



gsIges::gsIges(istream& is, gsXmlTree * data)
{
  m_xmldata = data;
  readIges(is);
}

//Function's return value is the directory entry pointer pushed back.
int MGIgesFstream::push_back_DE(MGIgesDirectoryEntry* de){
	int deNum=m_DirectoryEntries.size();
	m_DirectoryEntries.push_back(de);
	return deNum;
}
void gsIges::readIges(istream& is)
{
	m_StartSection=std::string();
	m_GSection=MGIgesGSec(filename);
	m_nlineGSec=0;
	m_DirectoryEntries.clear();
	MGIgesDirectoryEntry* de=new MGIgesDirectoryEntry;
	m_DirectoryEntries.push_back(de);//Set the dummy record.

	char lineData[73];
	char ID;
	int sequence;

	//Construct m_StartSection, m_GSection, m_DirectryEntry, m_ParamDataLines.

	//1. Start Section
	get_one_line(lineData,ID,sequence);
	if(!good())
		return;

	while(ID=='S'){
		m_StartSection+=lineData;
		get_one_line(lineData,ID,sequence);
		if(!good())
			return;
	}

	//2. Global Section
	std::string Gsecstring;
	m_nlineGSec=0;
	while(ID=='G'){
		Gsecstring+=lineData;
		m_nlineGSec++;
		get_one_line(lineData,ID,sequence);
		if(!good())
			return;
	}
	m_GSection.read_in(Gsecstring);

	//3. Directory Entry Section
	MGIgesDirectoryEntry* de;
	while(ID=='D'){
		std::string DEstring(lineData);
		get_one_line(lineData,ID,sequence); assert(ID=='D');
		DEstring+=lineData;//Concatenate two DE lines into one.
		de=new MGIgesDirectoryEntry(DEstring);
		m_DirectoryEntries.push_back(de);

		get_one_line(lineData,ID,sequence);
		if(!good())
			return;
	}

	char pDelimeter=m_GSection.paramDelimeter();//parameter delimeter

	//4. Parameter data Section
	while(ID=='P'){
		std::stringstream seq;seq.rdbuf()->str(lineData+65);
		//std::cout<<lineData<<std::endl;
		//std::cout<<seq.str()<<std::endl;
		int lnumber;seq>>lnumber;
		int DEpointer=MGIges::lnumber_to_DEpointer(lnumber);
		MGIgesDirectoryEntry& de=*(m_DirectoryEntries[DEpointer]);

		//Construct the string of the current Parameter data.
		int nlines=de.ParameterLineCount();
		std::string paramData(lineData,64);
		for(int i=1;i<nlines; i++){
			get_one_line(lineData,ID,sequence,64);
			paramData+=lineData;
			assert(ID=='P');
		}
		de.setPD(pDelimeter,paramData);

		//Get next line to process.
		get_one_line(lineData,ID,sequence);
		if(!good())
			return;
	}

	m_vertexListMap.set_ifstream(this);
	m_edgeListMap.set_ifstream(this);

//============================================================

    // Now, we need to scan the directory section
    // They are (as opposed to the global section) NOT stored in
    // a class variable for later use (corr. lines commented).
    // As we scan the D section we discover entities of type 128 and
    // type 126. They are immediately read into surf_ and crv_.
    // The P-numbers are also remembered.
    int sz = (int)sbufs[D].length();
    //cout << sz << endl << sbufs[D] << endl;
    int num_entries = sz/144;
    if ( num_entries*144 != sz)
      std::cerr<< "Directory section string has length != N*144, that is: a whole number times two lines.\n";
    // direntries_.resize(num_entries);
    const char* posD = sbufs[D].c_str();
    const char* posP0 = sbufs[P].c_str();
    //start_of_P_section_ = posP0;
    const char* posP = posP0;

    for (int i=0; i<num_entries; ++i) {

      IgesDirectory curr_entry(posD + i*144);

      // First we read all entities that may be included as part of
      // other entities (such as curve segments and surfaces, used for
      // composite curves and trimmed surfaces).
      // @@sbr We really should read all parts that are not created
      // using other entities.

      //direntries_[i] = curr_entry;

      // 	std::cout << i << ' ' << direntries_[i].entity_type_number << ' '
      // 	     << direntries_[i].param_data_start << ' '
      // 	     << direntries_[i].line_count << std::endl;
      // 	std::cout << (posP0 + 72*(direntries_[i].param_data_start-1)) << std::en//dl;

      posP = posP0 + 64*(curr_entry.param_data_start-1);

	int entity_number = curr_entry.entity_type_number;
	if (!supportedEntity(entity_number))
	{
	  std::cout<<"Unknown entity-type (" << entity_number <<
		    ") in file! Object neglected.\n";
	}
	else if (entity_number == 128)
        {
	  readIgesSurface(posP);
          //local_geom_.push_back(readIgesSurface(posP, curr_entry.line_count));
	  //local_colour_.push_back(curr_entry.color);
          //geom_id_.push_back(Pnumber_[i]);
          //geom_used_.push_back(0);
        }
	else if (entity_number == 126)
        {
	  readIgesCurve(posP);
          //local_geom_.push_back(crv);
	  //local_colour_.push_back(direntries_[i].color);
          //geom_id_.push_back(Pnumber_[i]);
          //geom_used_.push_back(0);
	  //pnumber_to_plane_normal_index_[Pnumber_[i]] = (int)plane_normal_.size()-1;
        }
	else if (entity_number == 116)
        {
	  readIgesPointCloud(posP);
          //local_geom_.push_back(pt_cl);
	  //local_colour_.push_back(direntries_[i].color);
          //geom_id_.push_back(Pnumber_[i]);
          //geom_used_.push_back(0);
        }
	else if (entity_number == 124)
        {
	  //std::tr1::shared_ptr< CoordinateSystem<3> > cs
	  // = readIGEStransformation(posP, direntries_[i].line_count);
	  //  coordsystems_[Pnumber_[i]] = *cs;
	}
    }
}

// Reads all objects of the IGES stored in the file
// as MGCL objects.
MGIgesIfstream& MGIgesIfstream::operator>>(MGGroup& group){
	int n=m_DirectoryEntries.size();//Num of directory entries.
	for(int i=1; i<n; i++){//i starts from 1 since 1st DE is dummy.
		MGIgesDirectoryEntry& de=*(m_DirectoryEntries[i]);
		if(!de.is_independent())
			continue;
		int tnum=de.EntityTypeNumber();
		if(tnum==TRANSFORMATION_MATRIX)//When transformation matrix
			continue;
		if(tnum==COLOR_DEFINITION)//When color definition
			continue;

		MGGel* gel=convert_to_gel(de);
		if(gel){
			//std::cout<<(*gel)<<std::endl;
			group.push_back(gel);
		}
	}
	return *this;
}


MGGel* MGIgesIfstream::convert_to_gel(
	const MGIgesDirectoryEntry& de
)const{
	int typeNumber=de.EntityTypeNumber();
	MGGel* gel=0;
	switch(typeNumber){
	case CIRCULAR_ARC: gel=convert_arc(de); break;
	case COMPOSITE_CURVE: gel=convert_composite(de); break;
	case CONIC_ARC: gel=convert_conic_arc(de); break;
	case PLANE: gel=convert_plane(de); break;
	case LINE: gel=convert_line(de); break;
	case PARAMETRIC_SPLINE_CURVE: gel=convert_spline(de); break;
	case MGIges::POINT: gel=convert_point(de); break;
	case RULED_SURFACE: gel=convert_ruled_surface(de); break;
	case SURFACE_OF_REVOLUTION: gel=convert_revolution_surface(de); break;
	case TABULATED_CYLINDER: gel=convert_tab_cyl(de); break;
	case RATIONAL_BSPLINE_CURVE: gel=convert_nurbs(de);break;
	case RATIONAL_BSPLINE_SURFACE: gel=convert_nurbs_surface(de); break;
	case BOUNDED_SURFACE: gel=convert_bounded_surface(de); break;
	case TRIMMED_SURFACE: gel=convert_trimmed_surface(de); break;
	case SPHERE: gel=convert_sphere158(de); break;
	case MANIFOLD_SOLID_BREP_OBJECT: gel=convert_MSBO(de); break;
	case PLANE_SURFACE: gel=convert_planeSurface(de); break;
	case RIGHT_CIRCULAR_CYLINDRICAL_SURFACE: gel=convert_cylinder(de); break;
	case SPHERICAL_SURFACE: gel=convert_sphere(de); break;
	case ASSOCIATIVITY_INSTANCE: gel=convert_group(de); break;
	case FACE: gel=convert_face(de); break;
	case MGIges::SHELL: gel=convert_shell(de); break;
	default:std::cout<<"MGIgesIfstream::convert_to_gel:Non object typeNumber:"
				<<typeNumber<<std::endl;
	}

	if(!gel)
		return 0;

	transform(de,*gel);
	MGAttribedGel* agel=dynamic_cast<MGAttribedGel*>(gel);

	//Visibility.
	if(!de.is_visible()){
		agel->set_no_display();
	}

	//Color
	int color=de.ColorNumber();
	MGColor* mcolor=0;
	if(color>0){
		mcolor=new MGColor(MGColor::get_instance(static_cast<MGColor::ColorID>(color)));
	}else if(color<0){
		const MGIgesDirectoryEntry* color_de=directoryEntry(-color);
		mcolor=convert_color(*color_de);
	}
	if(mcolor)
		agel->set_GLattrib(mcolor);

	//Line width
	int lw=de.LineWeightNumber();
	if(lw)
		agel->set_GLattrib(new MGLineWidth(de.LineWidth(GSection())));

	//Line Font
	int lf=de.LineFontPattern();
	if(lf>1){
		agel->set_GLattrib(new MGLineStipple(MGLineStipple::LineFont(lf)));
	}

	//std::cout<<(*gel)<<std::endl;////***********::
	return gel;
}

//Transform obj if de has the transformation matrix.
void MGIgesIfstream::transform(
	const MGIgesDirectoryEntry& de,	//de of the object obj.
	MGGel& obj					//Object to transform.
)const{
	int tid=de.transformID();
	if(!tid)
		return;

	const MGIgesDirectoryEntry& trde=*(m_DirectoryEntries[tid]);
	const MGIgesPD124* pd124=static_cast<const MGIgesPD124*>(trde.paramData().get());
	MGTransf tr;
	pd124->convert_to_MGTransf(tr);
	obj.transform(tr);
}


//From the current stream position, get one line data.
void MGIgesIfstream::get_one_line(
	char* lineData,	//line data without ID letter and sequence number
					//(that is, from column 1 to nchar) will be output.
					//buffer length must be >=(nchar+1).
	char& sectionID_letter,	//section identification letter of the line.
	int& sequence,	//ascending sequence number of the line.
	int nchar		//number of characters of one line
		//(When Parameter Data section nchar=64, and otherwise nchar=72)
){
	assert(nchar<=72);
	m_ifstream.get(lineData,nchar+1);
	char seqID[18];
	m_ifstream.get(seqID,74-nchar);
	sectionID_letter=seqID[72-nchar];
	m_ifstream>>sequence;//read sequence.
	char linefeed;
	m_ifstream.get(linefeed);	//read line feed.
}


bool gsIges::readSingleIgesLine(istream& is, char line_terminated[81],
				     IgesSection& sect)
//-----------------------------------------------------------------------------
{
    // Read any lonely endlines
    char c;
    while(is.get(c)) {
	if (c != '\n') {
	    is.putback(c);
	    break;
	}
    }

    // If we have reached end of file, return false
    if (is.eof()) return false;

    // We set the section indicator character to '\000' so
    // our switch further down is guaranteed to work.
    // Usually, this is unneeded, because the same is
    // done after the switch, but doing it here, too, 
    // means that we no longer depend on having a proper
    // line at the top of the file, and we are no longer
    // required to repeatedly call this function with the
    // same buffer as argument.
    line_terminated[72] = 0;

    // Read a line
    is.get(line_terminated, 81);

    // Get rid of the trailing newline, which would otherwise
    // terminate next get operation
    is.get(c);

    switch (line_terminated[72])
	{
	case 'S':
	    {
		sect = S;
		break;
	    }
	case 'G':
	    {
		sect = G;
		break;
	    }
	case 'D':
	    {
		sect = D;
		break;
	    }
	case 'P':
	    {
		sect = P;
		break;
	    }
	case 'T':
	    {
		sect = T;
		break;
	    }
	case '\000':
	  {
	    sect = E;
	    break;
	  }
	default:
            std::cerr<<"No valid section code for line.\n";
	}

    // Having read the section, we put a terminator there, so that
    // the line returned will be a string terminating just after the
    // content (including columns 0..71):
    line_terminated[72] = 0;

    return true;
}


void gsIges::readIgesCurve(const char* start)
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '126,.....' warn if not
    int type = readIGESint(start, pd, rd);
    assert (type==126);
    skipDelimiter(start, pd);

    // Read parameters
    int n1 = readIGESint(start, pd, rd) + 1;
    skipDelimiter(start, pd);
    int k1 = readIGESint(start, pd, rd) + 1;
    skipDelimiter(start, pd);
    int planar = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);
    int clo1;
    clo1 = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);
    int polynomial = readIGESint(start, pd, rd);
    int rational = (polynomial == 0);
    skipDelimiter(start, pd);
    int per1;
    per1 = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);

    int N = n1+k1;
    int i=0, j=0;
    gsVector<double> knot1(N);
    //cout << "\nknot1: \n";
    for (i=0; i<N; ++i) {
	knot1[i] = readIgesDouble(start, pd, rd);
	skipDelimiter(start, pd);
	//cout << knot1[i] << endl;
    }

    bool all_weights_are_one = true;
    N = n1;
    gsVector<double> weights(N);
    for (i=0; i<N; ++i) {
	weights[i] = readIgesDouble(start, pd, rd);
	if (weights[i] != 1.0) all_weights_are_one = false;
	skipDelimiter(start, pd);
    }

    N = n1*3;
    gsVector<double> coefs(N);
    for (i=0; i<N; ++i) {
	coefs[i] = readIgesDouble(start, pd, rd);
	skipDelimiter(start, pd);
    }

    double u1 = readIgesDouble(start, pd, rd);
    skipDelimiter(start, pd);
    double u2 = readIgesDouble(start, pd, rd);

    int dim = 3;

    //double norm[3];
    if (planar) // Discard normal vector of plane 
    {
      // Normal vector of plane
      for (i=0; i<3; ++i) {
	skipDelimiter(start, pd);
        //norm[i] = readIgesDouble(start, pd, rd);
	readIgesDouble(start, pd, rd);
      }
    }

    skipOptionalTrailingArguments(start, pd, rd);


    // NEW gsXmlNode
    // add ... gsGeometry
    


    // Extract the coordinate system  // DISABLED FOR NOW
    // int csentry = direntries_[direntry_index].trans_matrix;
    // CoordinateSystem<3> cs; 
    // if (csentry != 0) 
    // { // If value of directory entry is 0, we should
    //   // use identity.
    // 	map< int, CoordinateSystem<3> >::iterator it
    // 	    = coordsystems_.find(csentry);
    // 	if (it == coordsystems_.end()) {
    // 	    std::cout<<"Could not find the referred coordinate system ("
    // 		    << csentry << ") in the file. Using identity.\n";
    // 	} else {
    // 	    cs = it->second;
    // 	    std::cout<<"Transformation matrix for spline curve object "
    // 		    "is missing!\n";
    // 	}
    // }

//     // Need to skip ut to 5 arguments, as we can have
//     // both a normal (an error, but common) and extra
//     // property pointers.
//     //    skipOptionalTrailingArguments(start, pd, rd, 5);
//     // @afr: Instead, we just ignore everything until rd (usually semicolon).
//     while ((*start) != rd) {
// 	++start;
//     }
//     ++start;

/*
    SplineCurve *cvptr;
    if (rational == 1 && !all_weights_are_one)
    {
      int dim1 = dim+1;
      vector<double> coefs2(dim1*n1);
      for (i=0; i<n1; ++i) {
        for (j=0; j<dim; j++)
          coefs2[dim1*i+j] = coefptr[dim*i+j]*weights[i];
        coefs2[dim1*i+3] = weights[i];
      }
      cvptr = new SplineCurve(n1, k1, knot1.begin(), coefs2.begin(),
                              dim, true);
    }
    else
      cvptr = new SplineCurve(n1, k1, knot1.begin(), coefptr, dim);
    std::tr1::shared_ptr<SplineCurve> crv = (std::tr1::shared_ptr<SplineCurve>)(cvptr);

    double fuzzypar = 1e-10;
    if ((fabs(u1-crv->startparam()) > fuzzypar) ||
	(fabs(crv->endparam() - u2) > fuzzypar))
    {
#ifndef NDEBUG
        std::cout<<"Extracting subcurve for entity 126!\n";
	if (u1 < crv->startparam() || u2 > crv->endparam())
	  std::cout<<"Parameter value(s) outside domain, moving them inside!\n";
#endif
// #if 0
//         MESSAGE("Should extract subcurve! Skipping for now ...");
// #endif
	if (u1 < crv->startparam())
	  u1 = crv->startparam();
	if (u2 > crv->endparam())
	  u2 = crv->endparam();
	std::tr1::shared_ptr<SplineCurve> sub_crv(crv->subCurve(u1, u2));
	crv = sub_crv;
    }
*/

}

void gsIges::readIgesSurface(const char* start)
{
//Read in parameter data from string stream data.
void MGIgesPD128::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_integer(pDelimeter,pdstream,m_upper_indexU);
	get_integer(pDelimeter,pdstream,m_upper_indexV);
	get_integer(pDelimeter,pdstream,m_degreeU);
	get_integer(pDelimeter,pdstream,m_degreeV);

	get_integer(pDelimeter,pdstream,m_closedU);
	get_integer(pDelimeter,pdstream,m_closedV);
	get_integer(pDelimeter,pdstream,m_non_rational);
	get_integer(pDelimeter,pdstream,m_periodicU);
	get_integer(pDelimeter,pdstream,m_periodicV);

	int orderU=m_degreeU+1, orderV=m_degreeV+1;
	int nBrepU=m_upper_indexU+1, nBrepV=m_upper_indexV+1;

	//Read knot vector.
	int nBrepUPorderU=nBrepU+orderU, i;
	m_knotsU.size_change(orderU,nBrepU);
	for(i=0; i<nBrepUPorderU; i++)
		get_real(pDelimeter,pdstream,m_knotsU[i]);

	int nBrepVPorderV=nBrepV+orderV, j;
	m_knotsV.size_change(orderV,nBrepV);
	for(j=0; j<nBrepVPorderV; j++)
		get_real(pDelimeter,pdstream,m_knotsV[j]);
	
	int nBrepUBYnBrepV=nBrepU*nBrepV;

	//Read weights.
	m_weights.resize(nBrepU,nBrepV,1);
	for(j=0; j<nBrepV; j++){
		for(i=0; i<nBrepU; i++)
			get_real(pDelimeter,pdstream,m_weights(i,j,0));
	}

	//Read control points..
	m_control_points.resize(nBrepU,nBrepV,3);
	for(j=0; j<nBrepV; j++){
		for(i=0; i<nBrepU; i++){
			get_real(pDelimeter,pdstream,m_control_points(i,j,0));
			get_real(pDelimeter,pdstream,m_control_points(i,j,1));
			get_real(pDelimeter,pdstream,m_control_points(i,j,2));
		}
	}

	get_real(pDelimeter,pdstream,m_start_paramU);
	get_real(pDelimeter,pdstream,m_end_paramU);
	get_real(pDelimeter,pdstream,m_start_paramV);
	get_real(pDelimeter,pdstream,m_end_paramV);
}




// =====================================================


    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '128,.....' warn if not
    int type = readIGESint(start, pd, rd);
    assert( type==128);
    skipDelimiter(start, pd);

    // Read parameters
    int n1 = readIGESint(start, pd, rd) + 1;
    skipDelimiter(start, pd);
    int n2 = readIGESint(start, pd, rd) + 1;
    skipDelimiter(start, pd);
    int k1 = readIGESint(start, pd, rd) + 1;
    skipDelimiter(start, pd);
    int k2 = readIGESint(start, pd, rd) + 1;
    skipDelimiter(start, pd);
    int clo1, clo2;
    clo1 = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);
    clo2 = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);
    int polynomial = readIGESint(start, pd, rd);
    int rational = (polynomial == 0);
    skipDelimiter(start, pd);
    int per1, per2;
    per1 = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);
    per2 = readIGESint(start, pd, rd);
    skipDelimiter(start, pd);

    int N = n1+k1;
    int i=0, j=0;
    gsVector<double> knot1(N);
    //cout << "\nknot1: \n";
    for (i=0; i<N; ++i) {
	knot1[i] = readIgesDouble(start, pd, rd);
	skipDelimiter(start, pd);
	//cout << knot1[i] << endl;
    }

    N = n2+k2;
    //cout << "\nknot2: \n";
    gsVector<double> knot2(N);
    for (i=0; i<N; ++i) {
	knot2[i] = readIgesDouble(start, pd, rd);
	skipDelimiter(start, pd);
	//cout << knot2[i] << endl;
    }

    bool all_weights_are_one = true;
    N = n1*n2;
    gsVector<double> weights(N);
    for (i=0; i<N; ++i) {
	weights[i] = readIgesDouble(start, pd, rd);
	if (weights[i] != 1.0) all_weights_are_one = false;
	skipDelimiter(start, pd);
    }

    N = n1*n2;//*3
    gsMatrix<double> coefs(N,3);
    for (i=0; i<N; ++i) 
    {
      coefs(i,0) = readIgesDouble(start, pd, rd);
      skipDelimiter(start, pd);
      coefs(i,1) = readIgesDouble(start, pd, rd);
      skipDelimiter(start, pd);
      coefs(i,2) = readIgesDouble(start, pd, rd);
      skipDelimiter(start, pd);
    }

    double u1, u2, v1, v2;
    u1 = readIgesDouble(start, pd, rd);
    skipDelimiter(start, pd);
    u2 = readIgesDouble(start, pd, rd);
    skipDelimiter(start, pd);
    v1 = readIgesDouble(start, pd, rd);
    skipDelimiter(start, pd);
    v2 = readIgesDouble(start, pd, rd);

    skipOptionalTrailingArguments(start, pd, rd);
  
    

/*
    SplineSurface *sfptr;
    if (rational == 1 && !all_weights_are_one)
    {
      N = n1*n2;
      vector<double> coefs2(4*N);
      for (i=0; i<N; ++i) {
        for (j=0; j<3; j++)
          coefs2[4*i+j] = coefs[3*i+j]*weights[i];
        coefs2[4*i+3] = weights[i];
      }
      sfptr = new SplineSurface( n1, n2, k1, k2, knot1.begin(),
                                   knot2.begin(),
                                   coefs2.begin(), 3, true);
    }
    else
      sfptr = new SplineSurface( n1, n2, k1, k2, knot1.begin(),
                                   knot2.begin(),
                                   coefs.begin(), 3);
    std::tr1::shared_ptr<SplineSurface> srf(sfptr);
*/
}


void gsIges::readIgesPointCloud(const char* start)
{
    char pd = header_.pardel;
    char rd = header_.recdel;

    // Verify that start points to '116,.....' warn if not
    int type = readIGESint(start, pd, rd);
    assert( type==116 );
    skipDelimiter(start, pd);

    // Read parameters
    //Point pt(3);
    double pt[3];

    pt[0] = readIgesDouble(start, pd, rd);
    skipDelimiter(start, pd);
    pt[1] = readIgesDouble(start, pd, rd);
    skipDelimiter(start, pd);
    pt[2] = readIgesDouble(start, pd, rd);

    skipOptionalTrailingArguments(start, pd, rd);
/*
    std::tr1::shared_ptr<PointCloud3D> pt_cl(new PointCloud3D(pt.begin(), 1));
*/
}

  gsIges::IgesHeader::IgesHeader(const std::string gsec_string)
{
	std::istringstream gsstream(gsec_string);
	std::string onep;

	//1. parameter delimeter.
	char pdelim;
	gsstream.read(&pdelim,1);
	if(pdelim != ','){
		assert(pdelim=='1');
		char H;
		gsstream.read(&H,1); assert(H=='H');
		gsstream.read(&m_delimeter_param,1);
		gsstream.read(&pdelim,1);assert(pdelim==m_delimeter_param);
	}

	//2. record delimeter.
	get_Hollerith_string(m_delimeter_param,gsstream,onep);
	if(onep.length()){
		assert(onep.length()==1);
		m_delimeter_record=onep[0];
	}else
		m_delimeter_record=';';

	//3. productID_sender
	get_Hollerith_string(m_delimeter_param,gsstream,m_productID_sender);

	//4. file_name
	get_Hollerith_string(m_delimeter_param,gsstream,m_file_name);

	//5. native_systemID
	get_Hollerith_string(m_delimeter_param,gsstream,m_native_systemID);

	//6. m_preprocessor_version
	get_Hollerith_string(m_delimeter_param,gsstream,m_preprocessor_version);

	//7. m_number_of_bits_of_integer
	get_integer(m_delimeter_param,gsstream,m_number_of_bits_of_integer);

	//8. m_magnitude_single_precision
	get_integer(m_delimeter_param,gsstream,m_magnitude_single_precision);

	//9. m_significance_single_precision
	get_integer(m_delimeter_param,gsstream,m_significance_single_precision);

	//10. m_magnitude_double_precision
	get_integer(m_delimeter_param,gsstream,m_magnitude_double_precision);

	//11. m_significance_double_precision
	get_integer(m_delimeter_param,gsstream,m_significance_double_precision);

	//12. m_productID_receiver
	if(!get_Hollerith_string(m_delimeter_param,gsstream,m_productID_receiver))
		m_productID_receiver=m_productID_sender;

	//13. m_model_space_scale
	if(!get_real(m_delimeter_param,gsstream,m_model_space_scale))
		m_model_space_scale=1.;

	//14. m_unit_flag
	if(!get_integer(m_delimeter_param,gsstream,m_unit_flag))
		m_unit_flag=1;

	//15. m_unit_name
	if(!get_Hollerith_string(m_delimeter_param,gsstream,m_unit_name))
		m_unit_name="INCH";

	//16. m_max_number_of_line_weight_gradations
	if(!get_integer(m_delimeter_param,gsstream,m_max_number_of_line_weight_gradations))
		m_max_number_of_line_weight_gradations=1;

	//17. m_width_of_max_line_weight
	get_real(m_delimeter_param,gsstream,m_width_of_max_line_weight);

	//18. m_DateTime_File_generation
	get_Hollerith_string(m_delimeter_param,gsstream,m_DateTime_File_generation);

	//19. m_min_resolution
	get_real(m_delimeter_param,gsstream,m_min_resolution);

	//20. m_max_coordinate_value
	get_real(m_delimeter_param,gsstream,m_max_coordinate_value);

	//21. m_author_name
	get_Hollerith_string(m_delimeter_param,gsstream,m_author_name);

	//22. m_author_organazation
	get_Hollerith_string(m_delimeter_param,gsstream,m_author_organazation);

	//23. m_version_flag
	get_integer(m_delimeter_param,gsstream,m_version_flag);
	if(m_version_flag<1)
		m_version_flag=3;
	else if(m_version_flag>11)
		m_version_flag=11;

	//24. m_drafting_standard_flag
	get_integer(m_delimeter_param,gsstream,m_drafting_standard_flag);

	//25. m_DateTime_Model_generation
	get_Hollerith_string(m_delimeter_param,gsstream,m_DateTime_Model_generation);

	//26. m_application_protocolID
	get_Hollerith_string(m_delimeter_param,gsstream,m_application_protocolID);
}


gsIges::IgesDirectory::IgesDirectory(const char* start)
{
    char buffer[9];

    strncpy(buffer, start + 0, 9);
    entity_type_number = atoi(buffer);

    strncpy(buffer, start + 8, 9);
    param_data_start = atoi(buffer);

    strncpy(buffer, start + 16, 9);
    structure = atoi(buffer);

    strncpy(buffer, start + 24, 9);
    line_font_pattern = atoi(buffer);

    strncpy(buffer, start + 32, 9);
    level = atoi(buffer);

    strncpy(buffer, start + 40, 9);
    view = atoi(buffer);

    strncpy(buffer, start + 48, 9);
    trans_matrix = atoi(buffer);

    strncpy(buffer, start + 56, 9);
    buffer[8] = '\0'; // As the next element is not a space (which
		      // atoi expects).
    label_display = atoi(buffer);

    strncpy(buffer, start + 64, 9);
    buffer[8] = '\0';
    status = buffer;

    strncpy(buffer, start + 80, 9);
    line_weight = atof(buffer);  // @ DANGER, possible D-notation trouble

    strncpy(buffer, start + 88, 9);
    color = atoi(buffer);

    strncpy(buffer, start + 96, 9);
    line_count = atoi(buffer);

    strncpy(buffer, start + 104, 9);
    form = atoi(buffer);

    strncpy(buffer, start + 128, 9);
    buffer[8] = '\0';
    entity_label = buffer;

    strncpy(buffer, start + 136, 9);
    entity_number = atoi(buffer);
}

int gsIges::readIGESint(ccp& start, char pd, char rd)
{
  while (isspace(*start)) ++start;
  char numbuf[32];
  int numdig = 0;
  for(int i=0; i<31; i++) {
    if (start[i] == pd || start[i] == rd) {
      numdig = i;
      break;
    } else {
      numbuf[i] = start[i];
    }
  }
  numbuf[numdig] = 0; // Terminate numbuf
  start += numdig;
  return atoi(numbuf);
}


double gsIges::readIgesDouble(ccp& start, char pd, char rd)
{
    while (isspace(*start)) ++start;
    char numbuf[32];
    int numdig = 0;
    for(int i=0; i<31; i++) {
	if (start[i] == pd || start[i] == rd) {
	    numdig = i;
	    break;
	} else {
	    numbuf[i] = start[i];
	    if (numbuf[i] == 'D') numbuf[i] = 'E'; // FP notation...
	}
    }
    start += numdig; // Next value.
    // We must also disregard any trailing white spaces in numbuf.
    int nmb_trailing_spaces = 0;
    while (isspace(start[-1-nmb_trailing_spaces]))
	++nmb_trailing_spaces;
    numbuf[numdig-nmb_trailing_spaces] = 0; // Terminate numbuf

    std::stringstream ss ( numbuf );
    double val = -1.0;
    ss>> val;

    return val;
}

std::string gsIges::readIgesString(ccp& start, char delim, char delim2 )
{
    while (isspace(*start)) ++start;
    char numbuf[8];
    int numdig = 0;
    if (start[0] == delim || start[0] == delim2) {
	return string();
    }
    for(int i=0; i<7; i++) {
	if (start[i] == 'H') {
	    numdig = i;
	    break;
	} else {
	    numbuf[i] = start[i];
	}
    }
    numbuf[numdig] = 0; // Terminate numbuf
    //        cout << numdig << " " << numbuf << endl;
    int numchars = atoi(numbuf);
    if (numchars < 1) {
      std::cout<< "Less than one character in string: " << numchars<<".\n";
	return string();
    }
    start += numdig + 1 + numchars;
    //    cout << "[[[    " << numchars << "    ]]]" << endl;
    return string(start-numchars, numchars);
}

void gsIges::skipDelimiter(ccp& whereami, char pd)
{
    //  cout << "*whereami = " << *whereami << endl
    // << "pd        = " << pd << endl;

    // Skip whitespace
    while (isspace(*whereami)) ++whereami;
    if ((*whereami) == pd)
	++whereami;
    else {
	//char buf[101];
	//cout << endl << pd << endl;
	//cout << "-------------------------------------" << endl;
	//strncpy(buf,whereami,100);
	//cout << buf << endl;
      std::cout<< "Parsing anomaly, could not locate separator (" << pd << ") on p-line ..\n";
    }
}

bool gsIges::checkDelimiter(ccp& whereami, char wanted_delimiter,
                    char alternative_delimiter)
{
    //  cout << "*whereami = " << *whereami << endl
    // << "pd        = " << pd << endl;

    // Skip whitespace
    while (isspace(*whereami)) ++whereami;
    if ((*whereami) == wanted_delimiter)
      {
	++whereami;
	return true;
      }
    else if ((*whereami) == alternative_delimiter)
      {
	++whereami;
	return false;
      }
    else 
      {
	//char buf[101];
	//cout << endl << pd << endl;
	//cout << "-------------------------------------" << endl;
	//strncpy(buf,whereami,100);
	//cout << buf << endl;

	std::cout<< "Parsing anomaly, could not locate separator (" << wanted_delimiter 
		 << " or " << alternative_delimiter<< ") on p-line ..\n";
	return false;
      }
}

void gsIges::skipOptionalTrailingArguments(ccp& whereami,
					   char pd, char rd,
					   int max_to_skip)
{
  for (int j = 0; j < max_to_skip; ++j) {
    bool ismore = checkDelimiter(whereami, pd, rd);
    if (ismore) {
      // Read past the extra pointers in the file
      int npointer = readIGESint(whereami, pd, rd);
      int pdummy;
      for (int i = 0; i < npointer; i++) {
	skipDelimiter(whereami, pd);
	pdummy = readIGESint(whereami, pd, rd);
      }
    } else {
      break;
    }
  }
}



}// namespace internal

}// namespace gismo
