#ifdef EV_EXPORTS  

#define EV_API _declspec(dllexport) //表明标有此宏定义的函数和类是dll文件的导出函数和类，是dll文件的对外接口

#else

#define EV_API _declspec(dllimport) //表明标有此宏定义的函数和类的定义在dll文件中

#endif


#include "..\..\include\tnm\TNM_utility.h"
//#include "..\..\DllProj\transitNet\head.h"
#include "My_Predicate.h"


class ELINK;
class ENODE;
class ENET;
class EDest;
class EORIGIN;
class EPATH;
class DistanceDist;
class MULCLASS;
struct SubRange;
struct SubPath;
typedef map<int, SubPath*, less<int> > GroupPathMap;
typedef map<int, SubPath*, less<int> >::iterator GroupPathMapIter;


using namespace std;

class EV_API DistanceDist
{
public:
	DistanceDist() { ; };

	typedef enum {
		uni,
		lognormal,
		triangle,
	} Dists;

	DistanceDist(floatType min, floatType max, Dists tname)
	{
		lambda_min = min;
		lambda_max = max;
		DistType = tname;

		if (GetDistType() == DistanceDist::lognormal)
		{
			mean = 0.5*(lambda_min + lambda_max);
			sigma = 0.6;
			mu = (log(mean) + sigma * sigma / 2);
			delta_p = 0.5 * erf((log(lambda_max) - mu) / (sqrt(2) * sigma)) - \
				0.5 * erf((log(lambda_min) - mu) / (sqrt(2) * sigma));
		}

		if (GetDistType() == DistanceDist::triangle)
		{
			mean = 0.5 * (lambda_min + lambda_max);
		}
	}
	Dists				GetDistType() { return DistType;};
	floatType			GetCDF(floatType lambda);
	floatType			GetPDF(floatType lambda);

	Dists				DistType;
	floatType			lambda_min;
	floatType			lambda_max;
	floatType			mean;
	floatType			mu;
	floatType			sigma;
	floatType			delta_p;
};

struct EV_API EITERELEM
{
public:
	EITERELEM() {};
	int         iter; //current iteration number
	double      convRGap; //current relative gap
	double      stepsize; //current step size (rho)
	float       time; //current cputimeprint
	int			innerIters;
	int			numberofpaths;
	int			numberofsubpaths;
	int			numberofbuepath;
	float       ofv;
	float       numoftrips;
};

struct EV_API PATHELEM //HYPERLINK ELEM
{
	PATHELEM()
	{
		cost = 0.0;
		via = NULL;
	}
	floatType				cost;
	ELINK* via;// this is created for built strategy info
};





struct EV_API BSPLable
{
public:
	BSPLable(ENODE* v1, double v2, double v3, BSPLable* v4, ELINK* v5)
	{
		source = v1;
		objvalue1 = v2;
		objvalue2 = v3;
		preLable = v4;
		prelink = v5;
	};
	ENODE* source;
	BSPLable* preLable;
	ELINK* prelink;
	bool   nondominance;
	double objvalue1;//used for a more general bi-objective case,TIME 
	double objvalue2;//LENGTH
};

class EV_API ENODE
{
public:
	ENODE()
	{
		id = -1;
		xCord = 0;
		yCord = 0;
		scanStatus = 0;
		buffer = NULL;
		PathElem = new PATHELEM;
	};
	inline int		   id_() { return id; } //return id of this node.

	void			   ResetBSPLables() //earse all labels for a node
	{

		for (vector<BSPLable*>::iterator pv = labelSet.begin(); pv != labelSet.end(); pv++)
		{
			if (NULL != *pv)
			{
				delete* pv;
				*pv = NULL;
			}
		}
		labelSet.clear();
	}

	int                id;         //Node id
	largeInt           xCord;      //The x coordinate
	largeInt           yCord;      //The y coordinate
	vector<ELINK*>	   forwStar;   //All outgoing link pointers
	vector<ELINK*>	   backStar;   //All incoming link pointers
	floatType		   *buffer;    //It is a temporary working space for outside use.  
	smallInt           scanStatus; //A swith variable.  Normally it is used in shortest path computation, also somewhere else.
	bool               dummy;      //A switch variable, tell if the node is a temporary dummy node. In some applications, dummy nodes might be created.
	PATHELEM			*PathElem;  //Store the optimal routing policy.
	vector<BSPLable*>	labelSet;
};

class EV_API ELINK
{
public:
	ELINK()
	{
		id = 0;
		head = NULL;         /*starting node of the link*/
		tail = NULL;         /*ending node of the link */
		capacity = 0.0;          /*link capacity*/
		volume = 0.0;          /*link volume*/
		preFlow = 0.0;
		length = 0.0;          /*link length*/
		ffs = 0.0;          /*free flow speed*/
		fft = 0.0;          /*free flow travel time*/
		cost = 0;
		toll = 0.0;		
		fdCost = 0.0;
		markStatus = 0;
		buffer = NULL;
		m_isThrough = false;
	};


	floatType			GetCost(floatType v);
	floatType			GetIntCost();

	void				UpdateCost();
	void				UpdateDercost();


	
	inline int		    id_() { return id; } //return id of this link.

	int					id;			/*link id*/
	ENODE				*head;		/*starting node of the link*/
	ENODE				*tail;		/*ending node of the link */
	floatType			capacity;  /*link capacity, vph*/
	floatType			volume;    /*link volume, vph*/
	floatType			preFlow;  /*preload flow*/
	floatType			length;    /*link length, m*/
	floatType			ffs;       /*free flow speed, mph*/
	floatType			fft;       /*free flow travel time, h*/
	floatType           toll;      //this should be interpreted as any monetoary cost associated with link
	floatType			vot;		//value of time
	floatType			cost;      /*this can be a general concept of link cost (a basis for shortest path computation*/
	floatType           fdCost;    //this is the derivative of link cost;
	floatType			*buffer;    //this data is added as a working space for outside using.
	int					m_isThrough;//this is to mark a source link when expand the network 
	int					markStatus; //a status variable for temporary usage.
	bool				dummy; //if the link is a temporary "dummy" link.
	floatType			alpha; //parameters to define speed-flow relationship.
	floatType			beta;  //power
};

class EV_API EPATH
{
public:

	EPATH(EDest* s)
	{
		dest = s;
		flow = 0.0;
		cost = 0;
		markStatus = 0;
		fdCost = 0.0;
		preFlow = 0.0;
		preCost = 0.0;
		toll = 0.0;
		length = 0.0;
		newpath = true;
		id = 0;
	};
	void				  UpdatePathCost();
	void				  PrintPath();
	bool				  ChenckTopology();
	void				  SetPathid();


	EDest* dest;
	int                   id;
	floatType             flow;      // path flow
	floatType             preFlow;
	floatType             preCost;
	floatType             cost;      // additive path travel time
	floatType			  toll;      //additive toll 
	floatType			  length;
	floatType* buffer;   //a working space for outside using. 
	short                 markStatus;
	vector<ELINK*>		  vLinks;      //maintain path trace by link pointers
	floatType             fdCost; //first derviative
	floatType			  low_boundary;
	floatType			  up_boundary;
	bool				  newpath;
};

struct EV_API SubPath
{
public:
	SubPath(EPATH* p, floatType c)
	{
		path = p;
		subpathflow = c;
	}
	EPATH* path;
	floatType subpathflow;
};


struct EV_API SubRange
{
public:
	SubRange(EDest* d)
	{
		lowerBound = POS_INF_FLOAT;
		upperBound = POS_INF_FLOAT;
		dmd = 0;
		dest = d;
		addpath = false;
	}

	SubRange(floatType l, floatType u, floatType d,EDest* dd)
	{
		lowerBound = l;
		upperBound = u;
		dmd = d;
		dest = dd;
		addpath = false;
	}

	bool	Findpathbyid(int pid)
	{
		GroupPathMapIter it = GroupPaths.find(pid);
		if (it != GroupPaths.end())	return true;
		else
		{
			return false;
		}
	}

	void				UpdateRangePathSet();

	floatType			lowerBound;
	floatType			upperBound;
	floatType			dmd;
	
	//vector<EPATH*>		pathSet;//store paths from dest-based univeral path set
	//vector<floatType>		pathsubflow;//store subpath flows for each path
	
	GroupPathMap		GroupPaths;

	int					minIx;
	floatType			maxPathGap;
	floatType			currentRelativeGap;
	floatType           totalCost;
	EDest*				dest;
	bool				addpath;
};

class EV_API EDest
{
public:
	EDest()
	{
		dest = NULL;
		origin = NULL;
		assDemand = 0;
		buffer = NULL;
		costDif = 1.0;
		mid = -1;
		pid = 0;
		e_rs = 0.0;
	}

	ENODE					*dest;    //node on which the destination rest
	EORIGIN					*origin;
	floatType				assDemand;//assignment O-D demand
	floatType				upperDemand;//upper bound demand
	floatType				e_rs;		//excess demand flow
	floatType				gamma;		// for elastic demand function



	floatType* buffer; //this data is added as a working space for outside using.
	vector<EPATH*>			pathSet;
	floatType				costDif;
	vector<floatType>		pathSetbound;
	int						minIx; // this store the position for the minimum cost path;
	int						mid; // this store the index for od demand vector
	floatType				maxPathGap;
	floatType				currentRelativeGap;
	floatType				BBStepsize;//add for non-additive TAP problem

	vector<MULCLASS*>		classVector; // added for multi-class problem
	floatType				cof; //coefficient of this class determined by the problem, e.g., RTSM, VOT
	vector<SubRange*>		rangeVector;
	int						pid; // this record a id to each new generated path
	DistanceDist*			dDist;




	vector<EPATH*>			epaths;
	vector<floatType>		eboundary;
	void					ReleaseEMemory();


	floatType				DemandFunc(floatType u_rs); // D(u_rs)
	floatType				InverseDemandFunc(floatType dmd); // D-1(q)
	floatType				IntInverseDemandFunc(floatType dmd);//intergral of D-1(q)
	floatType				PrimeInverseDemandFunc(floatType dmd); //-1 * D-1'(q)
	floatType				IntWrs(floatType ers); //intergral of er



	bool            Initialize(EORIGIN* origin, ENODE* dest, floatType assDemand);
	inline int		id_();
	void			UpdatePathSetCost();
	void			RangeMerge();


	void			RemovePath(EPATH* path)//  delete this path from destination path set
	{
		for (vector<EPATH*>::iterator it = pathSet.begin(); it != pathSet.end(); it++)
		{
			if ((*it) == path)
			{
				pathSet.erase(it);
				break;
			}
		}
	}
	EPATH*	InPathSet(EPATH* path)
	{
		for (vector<EPATH*>::iterator it = pathSet.begin(); it != pathSet.end(); it++)
		{
			if ((*it)->vLinks == path->vLinks) return (*it);
		}
		return NULL;
	}

};



class EV_API EORIGIN
{
public:
	EORIGIN()  /*constructor*/
	{
		org = NULL;
		numOfDest = 0;
		destVector = NULL;
		m_tdmd = 0.0;
	};

	EORIGIN(ENODE* node, int nd)
	{
		org = node;
		numOfDest = nd;
		if (numOfDest > 0)
		{
			destVector = new EDest * [numOfDest];
			for (int i = 0; i < numOfDest; i++)   destVector[i] = new EDest;
		}
		else   destVector = NULL;
		m_tdmd = 0.0;
	};

	bool SetDest(int id, ENODE* node, floatType demand)
	{
		if (id > numOfDest || id <= 0)
		{
			cout << "\n\tSetDest in Origin Object: dest index exceeds the rannge!"
				<< "\n\trequired index = " << id << endl;;
			return false;
		}
		return destVector[id - 1]->Initialize(this, node, demand);
	};

	void SetDestMid(int id, int mid)
	{
		destVector[id - 1]->mid = mid;
	}

	inline int			id_()
	{
		return org->id;
	};

	ENODE				* org;			//node pointer
	int					numOfDest;     //number of destinations
	EDest				** destVector;  //destinations
	double				m_tdmd;      //total demand;
};

class EV_API MULCLASS
{
public:
	MULCLASS(EORIGIN* o, EDest* d, floatType i, floatType dmd)
	{
		classORG = o;
		classDest = d;
		dest = new EDest();
		dest->assDemand = dmd;
		dest->cof = i;
		dest->dest = d->dest;
	}

	EORIGIN* classORG;
	EDest* classDest;
	EDest* dest;
};

class EV_API ENET
{
public:
	ENET(const std::string& name)
	{
		networkName = name;
		numOfNode = 0;
		numOfLink = 0;
		numOfOrg = 0;
		numOfOD = 0;
		numOfTrips = 0;
		numOfPath = 0;
		numOfSubPath = 0;
		curIter = 0;
		numOfClass = 1;
		RGapIndicator = 1.0;
	};

	typedef enum {
		fd_single_gp,    // single pathj constraint
		fd_single_fw,
		fd_multiple_gp,	  // fixed demand for multi-class problem 
		fd_multiple_fw,
		fd_continuous_gp,	 // elastic demand for continuous problem
		fd_continuous_fw,
	}AlgType;


	// Functions
	int								BuildTAPASNet();
	void							ReadPathInfo();
	void							SetAlgorithm(AlgType algorithm) { ALG = algorithm; };
	string							GetAlgorithmName();
	int								Solver(AlgType algorithm);


	void							ReportIter();
	//void							ReportFlow();
	void							RecordCurrentIter();
	//void							RecordFlow();

	ENODE*							CatchNodePtr(int id, bool safe = false); //search from nodeVector by its id
	ELINK*							CatchLINKPtr(int id, bool safe = false); //search from nodeVector by its id
	ELINK*							CatchLinkPtr(ENODE* tail, ENODE* head);
	EORIGIN*						CatchORGPtr(int id, bool safe = false);
	EDest*							CatchDESTPtr(int id, EORIGIN* org);


	void							UpdateLinkNum() { numOfLink = (int)linkVector.size(); }
	void							UpdateNodeNum() { numOfNode = (int)nodeVector.size(); }
	EORIGIN* CreateOrigin(int nid, int nd); //create a static origin


	inline  bool					ReachAccuracy(double g) { return g <= convCriterion; }
	inline  bool					ReachMaxIter() { return curIter >= maxMainIter; } //test if maximum iteration is attained
	inline  bool					ReachMaxTime(floatType t) { return t >= maxIterTime; } /* running time (in min) */

	void							SetConv(floatType e) { convCriterion = e; };
	void							SetMaxIter(int i) { maxMainIter = i; };
	void							SetMaxIterTime(floatType i) { maxIterTime = i; };

	void							LabelSettingBSPTree(EORIGIN* Org);
	EPATH*							InitializePathByLabel(EORIGIN* Org, EDest* DEST ,BSPLable* label);




	// Solver for single-class distance constrainted TAP
	void							SolveSingleDCTAP_GP(); 
	void							UpdateNetLinkCost();
	void							InitialODFlow();
	void							SetODLengthLimit(floatType lenLimit);
	void							SinComputeRG();
	void							UpdateGPSubRangeflow(SubRange* range);
	void							DestColumnDropping(EDest* dest);

	// Solver for single-class distance constrainted TAP
	//void							SolveSingleDCTAP_FW();


	// Solver for multi-class distance constrainted TAP
	
	void							SolveMultiDCTAP_GP();
	void							InitialMultiClass();
	void							SetNumofClass(int i) { numOfClass = i; };
	void							MulComputeRG();
	void							InitialMultiFlow();
	void							RangeColumnDropping(SubRange* range);
	 
	

	 
	// Solver for continuous-class distance constrainted TAP
	void							SolveContinuousDCTAP_GP();
	void							Set_distance_Dist(floatType min, floatType max, DistanceDist::Dists type);
	void							InitialConODFlow();
	void							ConComputeRG();
	void							ConComputeRGII();
	void							LabelCorrectingSP(EORIGIN* org);
	void							IniRandomDistByFacor(floatType low_ratio, floatType up_ratio, DistanceDist::Dists type);
	void							IniRandomDistByLowUpValue(floatType low, floatType up, DistanceDist::Dists type);
	void							InfinityDistance(floatType low, floatType up, DistanceDist::Dists type);

	// Solver for continuous-class distance constrainted TAP with frank-wolfe
	void						    SolveContinuousDCTAP_FW();
	void							AllOrNothing();
	void							ConComputeRG_FW();
	void							InitialConODFlow_FW();
	void							LineSearch();
	void							UpdateFlowSolution();
	int								AllocateLinkBuffer(int size);


	AlgType							ALG;
	std::string						networkName;
	int								numOfClass;
	int								numOfNode;
	int								numOfLink;
	int								numOfOrg;
	int								numOfOD;		// number of O-D pairs
	floatType						numOfTrips;		// number of O-D pairs
	int								numOfPath;// number of paths
	int								curIter;       //current iteration
	int								numOfSubPath;
	floatType						OFV;			//objective function value
	int								innerIters;
	floatType* buffer;
	floatType						RGapIndicator; //convergence indicator
	floatType						convCriterion; // convergence gap
	floatType						errorCriterion; // convergence gap
	int								maxMainIter;
	floatType						maxIterTime;
	floatType						m_innerConv;
	clock_t							m_startRunTime; //this is the time when SolveTAP just called.
	floatType						stepsize;


	vector<ENODE*>					nodeVector;
	vector<ELINK*>					linkVector;
	vector<EORIGIN*>				orgVector;
	vector<EPATH*>					pathVector;
	vector<EDest*>					odVector;

	vector<EITERELEM*>				itersRecord;


};

