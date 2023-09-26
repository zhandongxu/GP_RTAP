#include <tchar.h>
#include "ev.h"

template <typename Map>
bool key_compare(Map const& lhs, Map const& rhs) {

	auto pred = [](auto a, auto b)
	{ return a.first == b.first; };

	if (lhs.size() == rhs.size()) {
		return std::equal(lhs.begin(), lhs.end(), rhs.begin(), pred);
	}
	else {
		return false;
	}
}


//sort time first and then LENGTH
bool lexicogr_order(const BSPLable* a, const BSPLable* b)
{
	if (a->objvalue1 < b->objvalue1) return true; //time asc, dis desc
	else if (a->objvalue1 == b->objvalue1 && a->objvalue2 < b->objvalue2) return true;
	else return false;
};

bool lengthdecreasing_order(const BSPLable* a, const BSPLable* b)
{
	if (a->objvalue2 > b->objvalue2) return true;
	else if (a->objvalue2 == b->objvalue2 && a->objvalue1 < b->objvalue1) return true;
	else return false;
}

floatType	DistanceDist::GetPDF(floatType lambda)
{
	floatType result = 0.0;
	if (lambda<lambda_min || lambda>lambda_max)
	{
		cout << "The lambda:" << lambda << " for pdf is out of range!!" << endl;
		system("PAUSE");
	}
	switch (GetDistType())
	{
	case DistanceDist::uni:
		result = 1.0 / (lambda_max - lambda_min);
		break;
	case DistanceDist::lognormal:
		result = 1.0 / (lambda * delta_p * sigma * sqrt(2 * PI)) * exp(-pow(log(lambda) - mu, 2) / (2 * sigma * sigma));
		break;
	case DistanceDist::triangle:
		if (lambda <= mean)
		{
			result = (lambda - lambda_min) * 4 / pow((lambda_max - lambda_min), 2);
		}
		else
		{
			result = (lambda_max - lambda) * 4 / pow((lambda_max - lambda_min), 2);
		}
		break;
	}
	return result;
}

floatType	DistanceDist::GetCDF(floatType lambda)
{
	floatType result = 0.0;
	if (lambda<lambda_min || lambda>lambda_max)
	{
		cout << "The lambda:" << lambda << " for pdf is out of range!!" << endl;
		system("PAUSE");
	}
	switch (GetDistType())
	{
	case DistanceDist::uni:
		result = (lambda - lambda_min) / (lambda_max - lambda_min);
		break;
	case DistanceDist::lognormal:
		result = 1.0 / (2 * delta_p) * (erf((log(lambda) - mu) / (sigma * sqrt(2))) - erf((log(lambda_min) - mu) / (sigma * sqrt(2))));
		break;
	case DistanceDist::triangle:
		if (lambda <= mean)
		{
			floatType tmp = pow(lambda_min, 2) * 2 / pow((lambda_max - lambda_min), 2) - lambda_min * lambda_min * 4 / pow((lambda_max - lambda_min), 2);
			result = pow(lambda,2) * 2/ pow((lambda_max - lambda_min), 2) - lambda_min * lambda * 4 / pow((lambda_max - lambda_min), 2) - tmp;
		}
		else
		{
			floatType tmp = lambda_max * mean * 4 / pow((lambda_max - lambda_min), 2) - pow(mean, 2) * 2 / pow((lambda_max - lambda_min), 2);
			result = lambda_max * lambda * 4 / pow((lambda_max - lambda_min), 2) - pow(lambda, 2) * 2 / pow((lambda_max - lambda_min), 2) + 0.5 - tmp;
		}
		break;
	}
	return result;

}

void	ELINK::UpdateDercost()
{

	if (volume <= 0.0)	fdCost = 0;
	else				fdCost = (fabs(beta - 0.0) < 1e-6 ? 0 : fft * alpha * beta * pow(volume / capacity, beta - 1.0) / capacity);

}

void	ELINK::UpdateCost()
{
	if (volume < -1e-8)
	{
		cout << "link:" << id << " has negiative flow:" << volume << endl;
		system("PAUSE");
	}
	else if (volume > -1e-8 && volume <= 0) volume = 0;

	if (volume == 0.0)	cost = fft;
	else				cost = fft * (1.0 + alpha * pow(volume / capacity, beta));

}

floatType	ELINK::GetIntCost()
{

	if (volume == 0.0)	return 0;
	else				return fft * volume * (1.0 + alpha * pow(volume / capacity, beta) / (beta + 1.0));
}

floatType ELINK::GetCost(floatType v)
{
	if (v < -1e-8)
	{
		cout << "link:" << id << " has negiative flow:" << v << endl;
		system("PAUSE");
	}
	else if (v > -1e-8 && v <= 0) v = 0;

	if (v == 0.0)	    return fft;
	else				return fft * (1.0 + alpha * pow(v / capacity, beta));
}


ENODE* ENET::CatchNodePtr(int id, bool safe)
{
	if (safe) return nodeVector[id - 1];
	else
	{
		vector<ENODE*>::iterator pv;
		pv = find_if(nodeVector.begin(), nodeVector.end(), predP(&ENODE::id_, id));
		if (pv == nodeVector.end()) return NULL;
		else                       return *pv;
	}
}

ELINK* ENET::CatchLinkPtr(ENODE* tail, ENODE* head)
{
	for (vector<ELINK*>::iterator pl = tail->forwStar.begin(); pl != tail->forwStar.end(); pl++)
	{
		if ((*pl))
		{
			if ((*pl)->head == head)
				return *pl;
		}
	}
	return NULL;
}

EORIGIN* ENET::CatchORGPtr(int id, bool safe)
{
	if (safe) return orgVector[id - 1];
	else
	{
		vector<EORIGIN*>::iterator pv;
		pv = find_if(orgVector.begin(), orgVector.end(), predP(&EORIGIN::id_, id));
		if (pv == orgVector.end()) return NULL;
		else                       return *pv;
	}
}

EORIGIN* ENET::CreateOrigin(int nid, int nd)
{
	ENODE* node = CatchNodePtr(nid);
	if (node == NULL)
	{
		cout << "\n\tnode " << nid << " is not a valid node object. " << endl;
		return NULL;
	}
	else
	{
		if (nd < 0)
		{
			cout << "\n\tOrigin " << node->id << " contains none or negative destinations." << endl;
			return NULL;
		}
		EORIGIN* org = new EORIGIN(node, nd);
		if (org == NULL)
		{
			cout << "\n\tCannot allocate memory for new origin" << endl;
			return NULL;
		}
		orgVector.push_back(org);
		numOfOrg++;
		numOfOD += nd;
		return org;
	}
}

void	SubRange::UpdateRangePathSet()
{
	double maxCost = -1.0, minCost = POS_INF_FLOAT;
	totalCost = 0.0;
	//double currentTotalCost = 0.0;
	double sumflow = 0;
	for (GroupPathMapIter it = GroupPaths.begin(); it != GroupPaths.end(); it++)
	{
		SubPath* subpath = it->second;
		EPATH* path = subpath->path;
		path->UpdatePathCost();
		if (path->cost > maxCost)	maxCost = path->cost;
		if (path->cost < minCost)
		{
			minCost = path->cost;
			minIx = it->first;
		}
		totalCost += path->cost * subpath->subpathflow;
		sumflow += subpath->subpathflow;
	}
	maxPathGap = maxCost - minCost;
	currentRelativeGap = fabs(1 - minCost * dmd / totalCost);
	//if (fabs(sumflow - dmd) > 1e-8)
	//{
	//	cout << "not convervative:" << fabs(sumflow - dmd) << endl;
	//}
}

bool EDest::Initialize(EORIGIN* o, ENODE* d, floatType as)
{
	if (o == NULL) return false;
	if (d == NULL) return false;
	origin = o;
	dest = d;
	assDemand = as;
	return true;

}

void EDest::UpdatePathSetCost()
{
	for (int i = 0; i < pathSet.size(); i++)
	{
		EPATH* path = pathSet[i];
		path->UpdatePathCost();
	}
}

void	EPATH::UpdatePathCost()
{
	cost = 0;
	for (int i = 0; i < vLinks.size(); ++i)
	{
		cost += vLinks[i]->cost;
	}
}

void	EPATH::SetPathid()
{
	id = dest->pid;
	dest->pid++;
}


int		ENET::BuildTAPASNet()
{
	string netfilename = networkName + "_net.dat";
	ifstream netFile;
	if (!TNM_OpenInFile(netFile, netfilename))    return 1;
	cout << "\tReading " << netfilename << "..." << endl;
	int nn, nl, nz, nft;

	TNM_SkipString(netFile, 3);
	netFile >> nz;
	TNM_SkipString(netFile, 3);
	netFile >> nn;
	TNM_SkipString(netFile, 3);
	netFile >> nft;
	TNM_SkipString(netFile, 3);
	netFile >> nl;
	string line;

	// we first create node objects
	for (int i = 0; i < nn; i++)
	{
		ENODE* node = new ENODE;
		node->id = i + 1;
		node->dummy = false;
		nodeVector.push_back(node);
	}

	int tail, head;
	floatType cap, len, fft, t, ffs, B, P, spd, toll;
	vector<string> words;
	int lid = 0;

	while (!netFile.eof())//for(int i =0 ;i<nl;i++)
	{
		getline(netFile, line);
		TNM_GetWordsFromLine(line, words);
		if (words.size() >= 1)
		{
			if (words[0].find_first_of("~") == -1 && words.size() >= 10)
			{
				lid++;
				TNM_FromString<int>(tail, words[0], std::dec);
				TNM_FromString<int>(head, words[1], std::dec);
				TNM_FromString<floatType>(cap, words[2], std::dec);
				TNM_FromString<floatType>(len, words[3], std::dec);
				TNM_FromString<floatType>(fft, words[4], std::dec);
				TNM_FromString<floatType>(B, words[5], std::dec);
				TNM_FromString<floatType>(P, words[6], std::dec);
				TNM_FromString<floatType>(spd, words[7], std::dec);
				TNM_FromString<floatType>(toll, words[8], std::dec);
				if (tail != head) //if tail = head, the link is ignored.
				{
					ENODE* tailnode = CatchNodePtr(tail, true);
					ENODE* headnode = CatchNodePtr(head, true);
					if (!tailnode || !headnode)
					{
						cout << "can not find tail and head of a link" << endl;
						return 2;
					}
					ELINK* link = new ELINK;
					link->tail = tailnode;
					link->head = headnode;
					link->id = lid;
					link->length = round(len*10)/10;//one decimal
					if (fft > 0)
					{			
						link->fft = fft; //in MIN
					}
					else
					{
						//link->length = 1.0;
						//link->ffs = 25.0;
						//link->fft = 0.001;
					}
					//check if capacity is abnormal, often centriod connectors's capacity are set to 0.0. 
					if (cap <= 1e-6)
					{
						cout << "Warning: link " << tail << " - " << head << "'s capacity = " << cap << ", its link performance function is forced to be constant." << endl;
						//B = 0.15;
						//P = 4;
					}
					if (toll < 0)
					{
						link->toll = 0.0;
					}
					else
					{
						link->toll = toll;
					}

					link->alpha = B;
					link->beta = P;
					link->capacity = cap;

					tailnode->forwStar.push_back(link);
					headnode->backStar.push_back(link);
					linkVector.push_back(link);

				}
				else
				{
					cout << "tail and head of a link are the same, check" << endl;
					system("PAUSE");
				}
			}
		}
	}

	UpdateLinkNum();
	UpdateNodeNum();
	// create O-D pair and trips
	int ix = 0;
	string odfilename = networkName + "_trp.dat";
	ifstream  odFile;
	if (!TNM_OpenInFile(odFile, odfilename))
	{
		return 2;
	}

	if (nz > 0)
	{
		cout << "\tReading " << odfilename << ", please wait..." << endl;
		for (int i = 0; i < 3; i++) getline(odFile, line);//this is to pass potential comments lines.
		string tStr;
		odFile >> tStr;

		while (tStr.compare("Origin") != 0 && !odFile.eof())
		{
			odFile >> tStr;
		}
		int orgID, nd, destID;
		floatType dmd;
		ENODE* node;
		vector<int> dvec;
		vector<floatType> dmdvec;
		EORIGIN* pOrg;
		for (int i = 0; i < nz; i++)
		{
			odFile >> orgID;
			if ((i + 1) % 100 == 0)
			{
				cout << "\t" << setw(4) << 100 * i / nz << "% completed" << endl;;
			}
			odFile >> tStr;
			nd = 0;
			if (!dvec.empty()) dvec.clear();
			if (!dmdvec.empty()) dmdvec.clear();
			while (tStr.compare("Origin") != 0 && !odFile.eof()) //not
			{
				if (TNM_FromString<int>(destID, tStr, std::dec))
				{
					odFile >> tStr;
					odFile >> tStr;
					TNM_FromString<floatType>(dmd, tStr, std::dec);
					odFile >> tStr;
					if (tStr == ";") odFile >> tStr; // test if ; has an space before it.
					//cout<<tStr<<endl;
					if (dmd > 0 && destID != orgID)
					{
						dvec.push_back(destID);
						dmdvec.push_back(dmd);
						nd++;
					}

				}
				else
				{
					cout << "destID = " << destID << endl;
					cout << "OD trip file includes unrecognized format!" << endl;
					return 5;
				}
			}
			if (dvec.size() > 0)
			{
				if ((pOrg = CreateOrigin(orgID, nd)) == NULL)
				{
					cout << "cannot create static origin object!" << endl;
					return 6;
				}
				pOrg->m_tdmd = 0.0;
				for (int j = 0; j < nd; j++)
				{
					node = CatchNodePtr(dvec[j], true);
					pOrg->SetDest(j + 1, node, dmdvec[j]);
					pOrg->SetDestMid(j + 1, ix); //set a postion number to a O-D pair
					ix++;
					pOrg->m_tdmd += dmdvec[j];
					numOfTrips += dmdvec[j];
				}
			}
		}
	}


	string tolfilename = networkName + "_tol.dat";
	ifstream tolfile;
	int	tolllinks = 0;
	if (TNM_OpenInFile(tolfile, tolfilename))
	{
		cout << "\tReading toll information from " << tolfilename << endl;
		for (int i = 0; i < numOfLink; i++) linkVector[i]->toll = 0.0;
		int ntl, ncount;
		floatType lowb, uppb;
		tolfile >> ntl >> lowb >> uppb;
		for (int i = 0; i < ntl; i++)
		{
			int from, to;
			floatType toll;
			tolfile >> from >> to >> toll;
			ENODE* fromnode, * tonode;
			if (from > 0 && from <= numOfNode && to > 0 && to <= numOfNode)
			{

				fromnode = CatchNodePtr(from, true);
				tonode = CatchNodePtr(to, true);
				if (from && to)
				{

					ELINK* link = CatchLinkPtr(fromnode, tonode);
					if (link)
					{
						if (toll < lowb)
						{
							link->toll = lowb;
						}
						else if (toll > uppb)
						{
							link->toll = uppb;
						}
						else
						{
							link->toll = toll;
							//ncount++;
						}
						if (link->toll > 0)	tolllinks++;

					}
					else
					{
						cout << "\tCannot locate link " << from << " - " << to << endl;

					}
				}
			}
			else
			{
				cout << "\tCannot locate link " << from << " - " << to << endl;
			}

			tolfile.close();
		}
	}

	cout << "===== node size:" << numOfNode << ",link size:" << numOfLink << ",od-pair size:" << numOfOD << ", trips:" << numOfTrips << ",numberoforg:" << numOfOrg << ",tolllinks:" << tolllinks << "=============" << endl;

	return 0;



}

string	ENET::GetAlgorithmName()
{
	string alg = "";
	switch (ALG)
	{
	case AlgType::fd_single_gp:
		alg = "fd_sin";
		break;
	case AlgType::fd_multiple_gp:
		alg = "fd_gp_multi" + to_string(numOfClass);
		break;
	case AlgType::fd_continuous_gp:
		alg = "fd_gp_con";
		break;
	case AlgType::fd_continuous_fw:
		alg = "fd_fw_con";
		break;
	default:
		cout << "No alg defined in the solver, please check solve!" << endl;
		return "";
		break;
	}
	return alg;
}

int		ENET::Solver(AlgType alg)
{
	SetAlgorithm(alg);
	cout << "Start running alg:" << GetAlgorithmName() << endl;
	switch (alg)
	{
	case AlgType::fd_single_gp:
		SolveSingleDCTAP_GP();
		break;
	//case AlgType::fd_single_fw:
		//SolveSingleDCTAP_FW();
		//break;
	case AlgType::fd_continuous_gp:
		SolveContinuousDCTAP_GP();
		break;
	case AlgType::fd_continuous_fw:
		SolveContinuousDCTAP_FW();
		break;
	case AlgType::fd_multiple_gp:
		SolveMultiDCTAP_GP();
		break;
	default:
		cout << "No alg defined in the solver, please check solve!" << endl;
		return 1;
		break;

	}
	return 0;
}

void	ENET::LabelSettingBSPTree(EORIGIN* Org)
{
	vector<BSPLable*> S;
	ENODE* node, * ownernode, * orgnode = Org->org;
	ELINK* link;
	for (int i = 0; i < numOfNode; i++)
	{
		node = nodeVector[i];
		node->ResetBSPLables();
	}
	BSPLable* newLable, * orgLable = new BSPLable(orgnode, 0, 0, NULL, NULL);
	orgnode->labelSet.push_back(orgLable);
	S.push_back(orgLable);
	while (S.size() > 0)
	{
		// sort lexicographically order first
		sort(S.begin(), S.end(), lexicogr_order);
		BSPLable* Lable = S[0];
		S.erase(S.begin());//remove the first label
		ownernode = Lable->source;

		for (vector<ELINK*>::iterator pv = ownernode->forwStar.begin(); pv != ownernode->forwStar.end(); pv++)
		{
			link = *pv;
			node = link->head;

			newLable = new BSPLable(node, Lable->objvalue1 + link->cost, Lable->objvalue2 + link->length, Lable, link);

			bool nondominate = true;

			for (vector<BSPLable*>::iterator it = node->labelSet.begin(); it != node->labelSet.end(); it++)
			{
				BSPLable* clabel = *it;
				if (fabs(newLable->objvalue1 - clabel->objvalue1) < 1e-10) clabel->objvalue1 = newLable->objvalue1;

				if (fabs(newLable->objvalue2 - clabel->objvalue2) < 1e-10) clabel->objvalue2 = newLable->objvalue2;

				if ((newLable->objvalue1 - clabel->objvalue1 > 0 && newLable->objvalue2 - clabel->objvalue2 >= 0) ||
					(newLable->objvalue1 - clabel->objvalue1 >= 0 && newLable->objvalue2 - clabel->objvalue2 > 0)
					)
				{
					nondominate = false; // this indicates newLable is dominated by clable
					break;
				}
			}
			if (nondominate)//newLabel not dominated by any of nb.labels(), check dominated labels in current set
			{
				vector<BSPLable*>::iterator it = node->labelSet.begin();
				while (it != node->labelSet.end())
				{
					BSPLable* clabel = *it;
					if ((clabel->objvalue1 - newLable->objvalue1 > 0 && clabel->objvalue2 - newLable->objvalue2 >= 0) ||
						(clabel->objvalue1 - newLable->objvalue1 >= 0 && clabel->objvalue2 - newLable->objvalue2 > 0))
					{
						it = node->labelSet.erase(it); //remove this dominated lable
						vector<BSPLable*>::iterator fit = find(S.begin(), S.end(), clabel);
						if (fit != S.end())	S.erase(fit);
					}
					else it++;
				}
				node->labelSet.push_back(newLable);
				S.push_back(newLable);
			}
		}
	}


}

EPATH* ENET::InitializePathByLabel(EORIGIN* org, EDest* DEST, BSPLable* label)
{
	BSPLable* destlabel = label;
	EPATH* path = new EPATH(DEST);
	while (label->source != org->org)
	{
		if (label->prelink)
		{
			path->vLinks.insert(path->vLinks.begin(), label->prelink);//add link to the front of link vec
			label = label->preLable;
		}
		else
		{
			cout << "Cannot find prelink of current minlable" << endl;
			system("PAUSE");
		}
	}
	if (path->vLinks.size() == 0)
	{
		cout << "The link size can not be zero!!" << endl;
		system("PAUSE");
	}
	path->cost = destlabel->objvalue1;
	path->length = destlabel->objvalue2;
	path->SetPathid();
	return path;
}

void	ENET::SetODLengthLimit(floatType lenLimit)
{
	EORIGIN* org;
	EDest* dest;
	for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		for (int j = 0; j < org->numOfDest; j++)
		{
			dest = org->destVector[j];
			SubRange* range = new SubRange(lenLimit, lenLimit,dest->assDemand,dest);
			dest->rangeVector.push_back(range);
		}
	}
}

void	ENET::InitialODFlow()
{
	EORIGIN* org;
	EDest* dest;
	SubRange* range;
	EPATH* path;
	for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		LabelSettingBSPTree(org);
		for (int j = 0; j < org->numOfDest; j++)
		{
			dest = org->destVector[j];
			for (vector<SubRange*>::iterator it = dest->rangeVector.begin(); it != dest->rangeVector.end(); it++)
			{
				SubRange* range = *it;
				sort(dest->dest->labelSet.begin(), dest->dest->labelSet.end(), lexicogr_order);
				for (int n = 0; n < dest->dest->labelSet.size(); n++)
				{
					BSPLable* label = dest->dest->labelSet[n];
					if (label->objvalue2 <= range->upperBound)
					{
						EPATH* newpath = InitializePathByLabel(dest->origin, dest , label);
						path = dest->InPathSet(newpath);
						if (path == NULL)
						{							
							dest->pathSet.push_back(newpath);
							SubPath* subpath = new SubPath(newpath, range->dmd);	
							pair<GroupPathMapIter, bool> ret = range->GroupPaths.insert(pair<int, SubPath*>(newpath->id, subpath));					
							newpath->flow += range->dmd;
							for (int k = 0; k < newpath->vLinks.size(); k++)
							{
								//cout << k << endl;;
								newpath->vLinks[k]->volume += range->dmd;
								newpath->vLinks[k]->UpdateCost();
								newpath->vLinks[k]->UpdateDercost();
							}
							break;
						}
					}
				}

			}
		}
	}

}

void	ENET::SolveSingleDCTAP_GP()
{
	m_startRunTime = clock();
	UpdateNetLinkCost();
	InitialODFlow();
	EORIGIN* org;
	EDest* dest;
	SubRange* range;

	if (!(ReachAccuracy(RGapIndicator) || ReachMaxIter() || ReachMaxTime((clock() - m_startRunTime) * 1.0 / (CLOCKS_PER_SEC * 60.0))))
	{
		do
		{
			curIter++;

			SinComputeRG();

			//perform inner loop
			innerIters = 0;
			int maxiters = 100;	// can also be adjusted
			while (innerIters < maxiters)
			{
				innerIters++;
				int badloop = 0;
				double maxgap = -1;
				//float scale = 1.0;
				for (int i = 0; i < numOfOrg; i++)
				{
					org = orgVector[i];
					for (int j = 0; j < org->numOfDest; j++)
					{
						dest = org->destVector[j];
						for (vector<SubRange*>::iterator it = dest->rangeVector.begin(); it != dest->rangeVector.end(); it++)
						{
							SubRange* range = *it;
							range->UpdateRangePathSet();
							if (range->currentRelativeGap > maxgap) maxgap = range->currentRelativeGap;
							if (range->GroupPaths.size() > 1 && range->currentRelativeGap > RGapIndicator)// can be adjusted by setting a "scale" parameter
							{
								badloop++;
								UpdateGPSubRangeflow(range);
								RangeColumnDropping(range);
							}
						}
						DestColumnDropping(dest);
					}
				}
				//
				if (badloop < 3)
				{
					//cout << maxgap << endl;
					break;
				}

			}

			RecordCurrentIter();

			//linkVector;

			cout << "iter:" << curIter << ",\t Relativegap:" << RGapIndicator << ",\t ttdmd:" << numOfTrips << ",\t OFV:" << OFV << ", \t innerIters:" << innerIters << endl;

		} while (!(ReachAccuracy(RGapIndicator) || ReachMaxIter() || ReachMaxTime((clock() - m_startRunTime) * 1.0 / (CLOCKS_PER_SEC * 60.0))));
	}
/*	for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		LabelCorrectingSP(org);
		for (int j = 0; j < org->numOfDest; j++)
		{
			dest = org->destVector[j];
			for (vector<EPATH*>::iterator it = dest->pathSet.begin(); it != dest->pathSet.end(); it++)
			{
	 			for (int k = 0; k < size((*it)->vLinks); k++)
				{
					if ((*it)->vLinks[k]->id == 1248 and (*it)->flow > 10)
					{
						cout << "path distance:" << (*it)->length << "flow:" << (*it)->flow << "Ori:" << org->org->id << "Dest:" << dest->dest->id << endl;
					}
				}
			}
		}
	}
*/
	
	for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		for (int j = 0; j < org->numOfDest; j++)
		{
			dest = org->destVector[j];
			if (org->org->id == 91 and dest->dest->id == 135)
			{
				for (vector<EPATH*>::iterator it = dest->pathSet.begin(); it != dest->pathSet.end(); it++)
				{
					for (int k = 0; k < size((*it)->vLinks); k++)
					{
						cout << "Path ID" << (*it)->id << "path distance:" << (*it)->length << "flow:" << (*it)->flow << "path cost:" << (*it)->cost <<
							"Link ID" << (*it)->vLinks[k]->id << endl;
					}
				}
			}
		}
	}

}


void	ENET::SolveContinuousDCTAP_GP()
{
	m_startRunTime = clock();
	UpdateNetLinkCost();
	InitialConODFlow();
	EORIGIN* org;
	EDest* dest;
	SubRange* range;

	if (!(ReachAccuracy(RGapIndicator) || ReachMaxIter() || ReachMaxTime((clock() - m_startRunTime) * 1.0 / (CLOCKS_PER_SEC * 60.0))))
	{
		do
		{
			curIter++;

			ConComputeRGII();

			//perform inner loop
			innerIters = 0;
			int maxiters = 500;	// can also be adjusted
			while (innerIters < maxiters)
			{
				innerIters++;
				int badloop = 0;
				double maxgap = -1;
				//float scale = 1.0;
				for (int i = 0; i < numOfOrg; i++)
				{
					org = orgVector[i];
					for (int j = 0; j < org->numOfDest; j++)
					{
						dest = org->destVector[j];
						for (vector<SubRange*>::iterator it = dest->rangeVector.begin(); it != dest->rangeVector.end(); it++)
						{
							SubRange* range = *it;
							range->UpdateRangePathSet();
							if (range->currentRelativeGap > maxgap) maxgap = range->currentRelativeGap;
							if (range->GroupPaths.size()>1 && range->currentRelativeGap > 0.1* RGapIndicator)
								// can be adjusted by setting a "scale" parameter
							{
								badloop++;
								UpdateGPSubRangeflow(range);
								RangeColumnDropping(range);
							}
						}

						DestColumnDropping(dest);
					}
				}
				if (badloop < 1)
				{
					//cout << maxgap << endl;
					break;
				}

			}
			
			RecordCurrentIter();

			cout << "iter:" << curIter << ",\t Relativegap:" << RGapIndicator << ",\t OFV:" << OFV << ", \t innerIters:" << innerIters << endl;

		} while (!(ReachAccuracy(RGapIndicator) || ReachMaxIter() || ReachMaxTime((clock() - m_startRunTime) * 1.0 / (CLOCKS_PER_SEC * 60.0))));
	}
	// print info for nd network
	/*for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		for (int j = 0; j < org->numOfDest; j++)
		{
			dest = org->destVector[j];
			cout << "org:" << org->org->id << ",dest:" << dest->dest->id << endl;
			for (vector<SubRange*>::iterator it = dest->rangeVector.begin(); it != dest->rangeVector.end(); it++)
			{
				SubRange* range = *it;
				cout <<"["<< range->lowerBound << "," << range->upperBound << "]" << endl;
				for (GroupPathMapIter itt = range->GroupPaths.begin(); itt != range->GroupPaths.end(); itt++)
				{
					SubPath* subpath = itt->second;
					cout << "len:"<<subpath->path->length << "," << subpath->subpathflow << endl;		
				}
			}
			cout << endl;
		}
	}*/
/*  
	for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		LabelCorrectingSP(org);
		for (int j = 0; j < org->numOfDest; j++)
		{
			dest = org->destVector[j];
			for (vector<EPATH*>::iterator it = dest->pathSet.begin(); it != dest->pathSet.end(); it++)
			{
				for (int k = 0; k < size((*it)->vLinks); k++)
				{
					if ((*it)->vLinks[k]->id == 1248 and (*it)->flow > 10)
					{
						cout << "path distance:" << (*it)->length << "flow:" << (*it)->flow << "Ori:" << org->org->id << "Dest:" << dest->dest->id << endl;
					}
				}
			}
		}
	}
// Disaggregate for a certain OD pair.
	*/ 
	/*for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		LabelCorrectingSP(org);
		for (int j = 0; j < org->numOfDest; j++)
		{
			dest = org->destVector[j];
			if (org->org->id == 91 and dest->dest->id == 135)
			{
				for (vector<EPATH*>::iterator it = dest->pathSet.begin(); it != dest->pathSet.end(); it++)
				{
					for (int k = 0; k < size((*it)->vLinks); k++)
					{
						cout << "Path ID" << (*it)->id << "flow:" << (*it)->flow << "path distance:" << (*it)->length << "path cost:" << (*it)->cost <<
							"Link ID" << (*it)->vLinks[k]->id << endl;
					}
				}
		
				cout << "Shortest path:" << dest->dest->PathElem->cost << endl;
			}
		}
	}
	*/
	ReadPathInfo();// read path info from files
}

void	ENET::ReadPathInfo()
{
	string filename = networkName + "_path.dat";
	ifstream infile;
	if (!TNM_OpenInFile(infile, filename))    return;
	cout << "\tReading the path file..." << endl;
	string pline;
	vector<string> words;

	if (!getline(infile, pline))
	{
		cout << "Failed to read the route informaiton from the following line: " << pline << endl;
		return ;
	}//skip the first line
	int n = 0;
	while (getline(infile, pline))
	{
		if (!pline.empty())//skip an empty line
		{
			vector<string> words,iwords;
			TNM_GetWordsFromLine(pline, words, ',', '"');
		
			TNM_GetWordsFromLine(words[1], iwords, '-', '"');
			
				
			ENODE *t,*h;
			double toll = 0;
			floatType len=0, cost=0;
			for (int i = 0; i < iwords.size()-1; i++)
			{
				int tid,hid;
				TNM_FromString(tid, iwords[i], std::dec);
				TNM_FromString(hid, iwords[i+1], std::dec);
				t = CatchNodePtr(tid);
				h = CatchNodePtr(hid);
				ELINK* link = CatchLinkPtr(t,h);
				len += link->length;
				cost += link->cost;
			}
			cout << words[0] << "," << len << "," << cost << endl;
			n++;
		}
		
	}
	numOfPath = n;
	cout << "\tRead in " << n << " paths" << endl;
	infile.close();



}

void	ENET::SolveContinuousDCTAP_FW()
{
	m_startRunTime = clock();
	UpdateNetLinkCost();
	AllocateLinkBuffer(1);
	InitialConODFlow_FW();
	EORIGIN* org;
	EDest* dest;
	SubRange* range;


	if (!(ReachAccuracy(RGapIndicator) || ReachMaxIter() || ReachMaxTime((clock() - m_startRunTime) * 1.0 / (CLOCKS_PER_SEC * 60.0))))
	{
		do
		{
			curIter++;

			AllOrNothing();
			LineSearch();
			UpdateFlowSolution();
			UpdateNetLinkCost();

			ConComputeRG_FW();
			RecordCurrentIter();
			cout << "iter:" << curIter << ",\t Relativegap:" << RGapIndicator << ",\t OFV:" << OFV << ", \t InnerIters:" << innerIters << endl;

		} while (!(ReachAccuracy(RGapIndicator) || ReachMaxIter() || ReachMaxTime((clock() - m_startRunTime) * 1.0 / (CLOCKS_PER_SEC * 60.0))));
	}

}

void	ENET::AllOrNothing()
{
	EORIGIN* org;
	EDest* dest;
	DistanceDist* dDist;

	for (int li = 0; li < numOfLink; li++)
	{
		linkVector[li]->buffer[0] = 0;  
	}

	for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		LabelSettingBSPTree(org);
		for (int j = 0; j < org->numOfDest; j++)
		{
			dest = org->destVector[j];
			dDist = dest->dDist;
			floatType rightbound = dDist->lambda_max;
			sort(dest->dest->labelSet.begin(), dest->dest->labelSet.end(), lexicogr_order);
			//dest->rangeVector;
			for (int n = 0; n < dest->dest->labelSet.size(); n++)
			{
				BSPLable* label = dest->dest->labelSet[n];
				if (label->objvalue2 < rightbound)
				{
					floatType lowerbound = __max(label->objvalue2, dDist->lambda_min);
					floatType dmd = (dDist->GetCDF(rightbound) - dDist->GetCDF(lowerbound)) * dest->assDemand;

					//intial new path
					EPATH* newpath = InitializePathByLabel(dest->origin, dest, label);
					for (int k = 0; k < newpath->vLinks.size(); k++)
					{
						newpath->vLinks[k]->buffer[0] += dmd;
					}
					rightbound = lowerbound;
					delete newpath;
				}
				if (rightbound == dDist->lambda_min) break;
			}
		}
	}


}

void	ENET::UpdateFlowSolution()
{
	for (int i = 0; i < numOfLink; i++)
	{
		linkVector[i]->volume = stepsize * linkVector[i]->buffer[0] + (1.0 - stepsize) * linkVector[i]->volume;// assign flow to link
	}
}

void	ENET::LineSearch()
{
	floatType dz;
	floatType a = 0;
	floatType b = 1.0;
	dz = 1.0;
	int iter = 0;
	while ((iter < 25 && (b - a) >= 1e-10 ) || dz > 0) //  
	{

		stepsize = (a + b) / 2.0;
		iter++;
		
		dz = 0.0;
		for (int i = 0; i < numOfLink; i++)
		{
			ELINK* link = linkVector[i];
			floatType tflow = link->volume + stepsize * (link->buffer[0] - link->volume);
			dz += (link->buffer[0] - link->volume) * link->GetCost(tflow);
		}

		if (dz < 0.0)
			a = stepsize;
		else
			b = stepsize;
	}
}

void	ENET::SolveMultiDCTAP_GP()
{
	m_startRunTime = clock();
	UpdateNetLinkCost();
	InitialMultiClass(); // 
	InitialMultiFlow(); // 

	EORIGIN* org;
	EDest *dest,*mDest;
	SubRange* range;
	if (!(ReachAccuracy(RGapIndicator) || ReachMaxIter() || ReachMaxTime((clock() - m_startRunTime) * 1.0 / (CLOCKS_PER_SEC * 60.0))))
	{
		do
		{
			curIter++;

			MulComputeRG(); 

			//perform inner loop
			innerIters = 0;
			int maxiters = 1000;
			while (innerIters < maxiters)
			{
				innerIters++;
				int badloop = 0;
				for (int i = 0; i < numOfOrg; i++)
				{
					org = orgVector[i];
					for (int j = 0; j < org->numOfDest; j++)
					{
						dest = org->destVector[j];
						for (int c = 0; c < numOfClass; c++)
						{
							mDest = dest->classVector[c]->dest;
							for (vector<SubRange*>::iterator it = mDest->rangeVector.begin(); it != mDest->rangeVector.end(); it++)
							{
								SubRange* range = *it;
								range->UpdateRangePathSet();

								if (range->currentRelativeGap > RGapIndicator / 1000.0)
								{
									badloop++;
									UpdateGPSubRangeflow(range); // flow equivalation	
								}
								if (innerIters % 5 == 0) RangeColumnDropping(range);//column dropping
							}

						}
					}
				}
				if (badloop < 1) break;
			}
			RecordCurrentIter();
			cout << "iter:" << curIter << ",\t Relativegap:" << RGapIndicator << ",\t OFV:" << OFV << ", \t innerIters:" << innerIters << endl;

		} 
		while (!(ReachAccuracy(RGapIndicator) || ReachMaxIter() || ReachMaxTime((clock() - m_startRunTime) * 1.0 / (CLOCKS_PER_SEC * 60.0))));
	}
	//EITERELEM* it = itersRecord.back();
	//cout << "time:" << TNM_FloatFormat(it->time, 20, 18);	
}

int		ENET::AllocateLinkBuffer(int size)
{
	ELINK* link;
	if (size <= 0)
	{
		cout << "\tInvalid size of link buffer array" << endl;
		return 1;
	}

	for (int i = 0; i < numOfLink; i++)
	{
		link = linkVector[i];
		if (link->buffer) delete[] link->buffer;
		link->buffer = new floatType[size];
		if (link->buffer == NULL)
		{
			cout << "\tCannot allocate memory for link buffer!" << endl;
			return 1;
		}
		for (int j = 0; j < size; j++)
			link->buffer[j] = 0.0;
	}
	return 0;
}

void	ENET::InitialMultiFlow()
{
	EORIGIN* org;
	EDest *dest, *mDest;
	SubRange* range;
	EPATH* path;
	for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		LabelSettingBSPTree(org);
		for (int j = 0; j < org->numOfDest; j++)
		{
			dest = org->destVector[j];
			sort(dest->dest->labelSet.begin(), dest->dest->labelSet.end(), lexicogr_order);
			for (int c = 0; c < numOfClass; c++)
			{
				mDest = dest->classVector[c]->dest;
				for (vector<SubRange*>::iterator it = mDest->rangeVector.begin(); it != mDest->rangeVector.end(); it++)
				{
					SubRange* range = *it;
					for (int n = 0; n < dest->dest->labelSet.size(); n++)
					{
						BSPLable* label = dest->dest->labelSet[n];
						if (label->objvalue2 <= range->upperBound)
						{
							EPATH* newpath = InitializePathByLabel(dest->origin, dest, label);
							path = dest->InPathSet(newpath);
							if (path == NULL)
							{
								dest->pathSet.push_back(newpath);
								SubPath* subpath = new SubPath(newpath, range->dmd);
								pair<GroupPathMapIter, bool> ret = range->GroupPaths.insert(pair<int, SubPath*>(newpath->id, subpath));
								newpath->flow += range->dmd;
								for (int k = 0; k < newpath->vLinks.size(); k++)
								{
									newpath->vLinks[k]->volume += range->dmd;
									newpath->vLinks[k]->UpdateCost();
									newpath->vLinks[k]->UpdateDercost();
								}
								
							}
							else
							{
								SubPath* subpath = new SubPath(path, range->dmd);
								pair<GroupPathMapIter, bool> ret = range->GroupPaths.insert(pair<int, SubPath*>(path->id, subpath));		
								path->flow += range->dmd;
								for (int k = 0; k < path->vLinks.size(); k++)
								{
									path->vLinks[k]->volume += range->dmd;
									path->vLinks[k]->UpdateCost();
									path->vLinks[k]->UpdateDercost();
								}		
							}
							break;
						}
					}

				}
			}		
		}
	}
}


void	ENET::InitialMultiClass()
{
	for (int i = 0; i < numOfOrg; i++)
	{
		EORIGIN* org = orgVector[i];
		for (int j = 0; j < org->numOfDest; j++)
		{
			EDest* dest = org->destVector[j];
			//First obtain a set of distance
			floatType interval = (dest->dDist->lambda_max - dest->dDist->lambda_min) / (numOfClass * 1.0);
			vector<floatType> lowLambdas, Lambdas;
			floatType lambda;
			lowLambdas.push_back(dest->dDist->lambda_min);

			if (numOfClass == 1)
			{
				lambda = (dest->dDist->lambda_max + dest->dDist->lambda_min) / 2.0;
				Lambdas.push_back(lambda);
				lowLambdas.push_back(dest->dDist->lambda_max);
			}
			else
			{
				lambda = dest->dDist->lambda_min + interval / 2.0;
				Lambdas.push_back(lambda);
				lowLambdas.push_back(lambda + interval / 2.0);
				int ix = 2;
				while (ix <= numOfClass)
				{
					lambda =  Lambdas.back() + interval;
					lowLambdas.push_back(lambda + interval / 2.0);
					Lambdas.push_back(lambda);
					ix = ix + 1;
				}
			}

			if (Lambdas.size() != numOfClass)
			{
				cout << "The number of class should be equal to the size of betas, please check!" << endl;
				system("PAUSE");
			}

			lowLambdas[0] = dest->dDist->lambda_min;
			lowLambdas.back() = dest->dDist->lambda_max;

			for (int n = 0; n < numOfClass; n++)
			{
				floatType demand = dest->assDemand * (dest->dDist->GetCDF(lowLambdas[n + 1]) - dest->dDist->GetCDF(lowLambdas[n]));
				MULCLASS* cls = new MULCLASS(org, dest, Lambdas[n], demand);
				SubRange* range = new SubRange(Lambdas[n], Lambdas[n], demand, dest); // one value
				cls->dest->rangeVector.push_back(range);
				dest->classVector.push_back(cls);
			}
		}
	}
}


void	ENET::InitialConODFlow_FW()
{
	EORIGIN* org;
	EDest* dest;
	DistanceDist* dDist;
	for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		LabelSettingBSPTree(org);
		for (int j = 0; j < org->numOfDest; j++)
		{
			dest = org->destVector[j];
			dDist = dest->dDist;
			floatType rightbound = dDist->lambda_max;
			sort(dest->dest->labelSet.begin(), dest->dest->labelSet.end(), lexicogr_order);
			//cout << dDist->GetCDF(dDist->lambda_min) <<"，"<< dDist->GetCDF(dDist->lambda_max)<<endl;
			for (int n = 0; n < dest->dest->labelSet.size(); n++)
			{
				BSPLable* label = dest->dest->labelSet[n];
				if (label->objvalue2 < rightbound)
				{
					floatType lowerbound = __max(label->objvalue2, dDist->lambda_min);
					floatType dmd = (dDist->GetCDF(rightbound) - dDist->GetCDF(lowerbound)) * dest->assDemand;
					
					//intial new path
					EPATH* newpath = InitializePathByLabel(dest->origin, dest, label);			
					for (int k = 0; k < newpath->vLinks.size(); k++)
					{
						newpath->vLinks[k]->volume += dmd;
						newpath->vLinks[k]->UpdateCost();
						newpath->vLinks[k]->UpdateDercost();
					}
					rightbound = lowerbound;
					delete newpath;
				}
				if (rightbound == dDist->lambda_min) break;
			}
		}
	}
}


void	ENET::InitialConODFlow()
{
	EORIGIN* org;
	EDest* dest;
	//SubRange* range;
	//EPATH* path;
	DistanceDist* dDist;
	for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		LabelSettingBSPTree(org);
		for (int j = 0; j < org->numOfDest; j++)
		{
			dest = org->destVector[j];
			dDist = dest->dDist;
			floatType rightbound = dDist->lambda_max;
			sort(dest->dest->labelSet.begin(), dest->dest->labelSet.end(), lexicogr_order);
			//dest->rangeVector;
			for (int n = 0; n < dest->dest->labelSet.size(); n++)
			{
				BSPLable* label = dest->dest->labelSet[n];
				if (label->objvalue2 < rightbound)
				{
					floatType lowerbound = __max(label->objvalue2, dDist->lambda_min);
					floatType dmd = (dDist->GetCDF(rightbound) - dDist->GetCDF(lowerbound)) * dest->assDemand;
					SubRange* range = new SubRange(lowerbound, rightbound, dmd, dest);
					dest->rangeVector.push_back(range);

					//intial new path
					EPATH* newpath = InitializePathByLabel(dest->origin, dest, label);
					dest->pathSet.push_back(newpath);
					SubPath* subpath = new SubPath(newpath, range->dmd);
					pair<GroupPathMapIter, bool> ret = range->GroupPaths.insert(pair<int, SubPath*>(newpath->id, subpath));
					newpath->flow += range->dmd;
					for (int k = 0; k < newpath->vLinks.size(); k++)
					{
						newpath->vLinks[k]->volume += range->dmd;
						newpath->vLinks[k]->UpdateCost();
						newpath->vLinks[k]->UpdateDercost();
					}
					rightbound = lowerbound;
				}
				if (rightbound == dDist->lambda_min) break;
			}
			//reverse(dest->rangeVector.begin(), dest->rangeVector.end()); // in an increasing order
		}
	}
}

void	ENET::Set_distance_Dist(floatType min, floatType max, DistanceDist::Dists type)
{
	EORIGIN* org;
	EDest* dest;
	for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		for (int j = 0; j < org->numOfDest; j++)
		{
			dest = org->destVector[j];
			DistanceDist* dist = new DistanceDist(min, max, type);
			dest->dDist = dist;
		}
	}

}

void	ENET::UpdateGPSubRangeflow(SubRange* range)
{

	if (range->GroupPaths.size() <= 1) return;
	SubPath* minsubpath = range->GroupPaths[range->minIx];
	EPATH* MinCostpath = minsubpath->path;
	ELINK* link;
	floatType dflow;

	for (GroupPathMapIter it = range->GroupPaths.begin(); it != range->GroupPaths.end(); it++)
	{
		SubPath* subpath = it->second;
		EPATH* path = subpath->path;
		if (path != MinCostpath)
		{
			path->UpdatePathCost();
			MinCostpath->UpdatePathCost();

			floatType DerSum = 0.0;
			for (int i = 0; i < MinCostpath->vLinks.size(); i++)
			{
				link = MinCostpath->vLinks[i];
				DerSum += link->fdCost;
				link->markStatus = 1;// tag 1 for arcs on the shortest path
			}

			for (int i = 0; i < path->vLinks.size(); i++)
			{
				link = path->vLinks[i];
				if (link->markStatus == 1) DerSum -= link->fdCost;
				else
				{
					DerSum += link->fdCost;
				}
			}
			if (DerSum <= 0)
			{
				DerSum = 1e-10;
				//cout << "The second order derivatives is negative:" << DerSum << ", please check!" << endl;
			}

			double dev = path->cost - MinCostpath->cost;
			double dflow = 0.0;
			if (DerSum > 0.0)
			{
				if (dev >= 0)
				{
					dflow = __min(1.0 * dev / DerSum, subpath->subpathflow) * 0.95;
				}
				else
				{
					dflow = __max(1.0 * dev / DerSum, -minsubpath->subpathflow) * 0.95;
				}
			}
			subpath->subpathflow -= dflow;		
			subpath->path->flow -= dflow;
			minsubpath->subpathflow += dflow;
			minsubpath->path->flow += dflow;

			if (fabs(path->flow) < 1e-10) path->flow = 0;
			if (fabs(MinCostpath->flow) < 1e-10) MinCostpath->flow = 0;

			//if (path->flow < 0 || MinCostpath->flow < 0)
			//{
			//	cout << "negative path flow:"<< path->flow << "," << MinCostpath->flow << endl;
			//	system("PAUSE");
			//}

			for (int i = 0; i < MinCostpath->vLinks.size(); i++)
			{
				MinCostpath->vLinks[i]->volume += dflow;
				MinCostpath->vLinks[i]->UpdateCost();
				MinCostpath->vLinks[i]->UpdateDercost();
				MinCostpath->vLinks[i]->markStatus = 0;
			}
			for (int i = 0; i < path->vLinks.size(); i++)
			{
				path->vLinks[i]->volume -= dflow;
				path->vLinks[i]->UpdateCost();
				path->vLinks[i]->UpdateDercost();
				path->vLinks[i]->markStatus = 0;
			}
		}
		
	}

}

void	ENET::DestColumnDropping(EDest* dest)
{
	EPATH* maxpath = NULL;
	double maxflow = -1.0, ttflow = 0.0, dflow;
	for (vector<EPATH*>::iterator it = dest->pathSet.begin(); it != dest->pathSet.end();)
	{
		EPATH* path = *it;
		if (maxflow < path->flow)
		{
			maxflow = path->flow;
			maxpath = *it;
		}

		if (fabs(path->flow) < 1e-10)
		{
			for (vector<SubRange*>::iterator ir = dest->rangeVector.begin(); ir != dest->rangeVector.end(); ir++)
			{
				SubRange* range = *ir;
				GroupPathMapIter itt = range->GroupPaths.begin();
				while (itt != range->GroupPaths.end())
				{
					if (itt->second->path == path)
					{
						itt = range->GroupPaths.erase(itt);
					}
					else itt++;
				}
				/*int n = range->pathSet.size() - 1;
				while (n > 0)
				{
					EPATH* cpath = range->pathSet[n];
					if (cpath == path)
					{
						range->pathSet.erase(range->pathSet.begin() + n);
						range->pathsubflow.erase(range->pathsubflow.begin() + n);
					}
					n--;
				}*/
			}

			delete* it;
			it = dest->pathSet.erase(it);
		}
		else
		{
			ttflow += path->flow;
			it++;
		}
	}

	dflow = dest->assDemand - ttflow;

	//if (fabs(dflow) > 1e-8)
	//{
	//	cout << "Unconservation flow :" << fabs(dflow) << " for org:" << dest->origin->org->id << ", dest:" << dest->dest->id << ", be careful!" << endl;
	//}

}

void	ENET::RangeColumnDropping(SubRange* range)
{
	GroupPathMapIter itt = range->GroupPaths.begin();
	while (itt != range->GroupPaths.end())
	{
		if (fabs(itt->second->subpathflow) < 1e-10)
		{
			itt = range->GroupPaths.erase(itt);
		}
		else itt++;
	}
}

void	ENET::UpdateNetLinkCost()
{
	for (int i = 0; i < numOfLink; i++)
	{
		ELINK* link = linkVector[i];
		link->UpdateCost();
		link->UpdateDercost();
	}
}

void	ENET::SinComputeRG()
{
	EORIGIN* org;
	EDest* dest;
	SubRange* range;
	EPATH* path;
	OFV = 0;
	numOfPath = 0;
	floatType   ttcost = 0.0;
	floatType	ttmincost = 0.0;

	for (int i = 0; i < numOfLink; i++)
	{
		OFV += linkVector[i]->GetIntCost();
	}

	for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		LabelSettingBSPTree(org);
		for (int j = 0; j < org->numOfDest; j++)
		{
			dest = org->destVector[j];
			sort(dest->dest->labelSet.begin(), dest->dest->labelSet.end(), lexicogr_order);
			dest->UpdatePathSetCost();
			numOfPath += dest->pathSet.size();
			for (vector<SubRange*>::iterator it = dest->rangeVector.begin(); it != dest->rangeVector.end(); it++)
			{
				SubRange* range = *it;
				for (GroupPathMapIter it = range->GroupPaths.begin(); it != range->GroupPaths.end(); it++)
				{
					ttcost += it->second->subpathflow * it->second->path->cost;
				}
				
				for (int n = 0; n < dest->dest->labelSet.size(); n++)
				{
					BSPLable* label = dest->dest->labelSet[n];
					if (label->objvalue2 <= range->upperBound) //uppper bound: dist limit
					{
						ttmincost += label->objvalue1 * range->dmd;

						EPATH* newpath = InitializePathByLabel(dest->origin, dest, label);
						path = dest->InPathSet(newpath);
						if (path == NULL)
						{
							dest->pathSet.push_back(newpath);
							SubPath* subpath = new SubPath(newpath, 0.0);
							pair<GroupPathMapIter, bool> ret = range->GroupPaths.insert(pair<int, SubPath*>(newpath->id, subpath));
						}
						else
						{
							// check whether path has been included into the subrange path set
							if (!range->Findpathbyid(path->id))
							{
								SubPath* subpath = new SubPath(path, 0.0);
								pair<GroupPathMapIter, bool> ret = range->GroupPaths.insert(pair<int, SubPath*>(path->id, subpath));
							}
						}
						break;
					}
				}
			}
		}
	}

	RGapIndicator = fabs(1.0 - ttmincost / ttcost);
}

void	ENET::MulComputeRG()
{
	EORIGIN* org;
	EDest* dest, * mDest;
	SubRange* range;
	EPATH* path;
	OFV = 0;
	numOfPath = 0;
	numOfSubPath = 0;
	floatType   ttcost = 0.0;
	floatType	ttmincost = 0.0;
	for (int i = 0; i < numOfLink; i++)
	{
		OFV += linkVector[i]->GetIntCost();
	}
	
	for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		LabelSettingBSPTree(org);
		for (int j = 0; j < org->numOfDest; j++)
		{
			numOfOD += org->numOfDest;
			dest = org->destVector[j];
			sort(dest->dest->labelSet.begin(), dest->dest->labelSet.end(), lexicogr_order);
			dest->UpdatePathSetCost();
			numOfPath += dest->pathSet.size();
			int a = 0;
			for (int c = 0; c < numOfClass; c++)
			{
				mDest = dest->classVector[c]->dest;
				for (vector<SubRange*>::iterator it = mDest->rangeVector.begin(); it != mDest->rangeVector.end(); it++)
				{
					SubRange* range = *it;
					for (GroupPathMapIter it = range->GroupPaths.begin(); it != range->GroupPaths.end(); it++)
					{
						ttcost += it->second->subpathflow * it->second->path->cost;
					}
					for (int n = 0; n < dest->dest->labelSet.size(); n++)
					{
						
						BSPLable* label = dest->dest->labelSet[n];
						if (label->objvalue2 <= range->upperBound)
						{
							ttmincost += label->objvalue1 * range->dmd;
							EPATH* newpath = InitializePathByLabel(dest->origin, dest, label);
							path = dest->InPathSet(newpath);
							if (path == NULL)
							{
								dest->pathSet.push_back(newpath);
								SubPath* subpath = new SubPath(newpath, 0.0);
								pair<GroupPathMapIter, bool> ret = range->GroupPaths.insert(pair<int, SubPath*>(newpath->id, subpath));
							}
							else
							{
								// check whether path has been included into the subrange path set
								if (!range->Findpathbyid(path->id))
								{
									SubPath* subpath = new SubPath(path, 0.0);
									pair<GroupPathMapIter, bool> ret = range->GroupPaths.insert(pair<int, SubPath*>(path->id, subpath));
								}
							}
							break;
						}
					}	
				}
				numOfSubPath += mDest->rangeVector.size();
			}
		}
	}
	cout << "subNum:" << numOfSubPath << endl;
	RGapIndicator = fabs(1.0 - ttmincost / ttcost);
}

void	ENET::ConComputeRG_FW()
{
	EORIGIN* org;
	EDest* dest;
	DistanceDist* dDist;
	OFV = 0;
	numOfPath = 0;
	floatType   ttcost = 0.0;
	floatType	ttmincost = 0.0;
	EPATH* path;

	for (int i = 0; i < numOfLink; i++)
	{
		OFV += linkVector[i]->GetIntCost();
		ttcost += linkVector[i]->volume * linkVector[i]->cost;
	}

	for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		LabelSettingBSPTree(org);
		for (int j = 0; j < org->numOfDest; j++)
		{
			dest = org->destVector[j];
			dDist = dest->dDist;
			floatType rightbound = dDist->lambda_max;
			bool hitlowerbound = false;
			sort(dest->dest->labelSet.begin(), dest->dest->labelSet.end(), lexicogr_order);

			for (int n = 0; n < dest->dest->labelSet.size(); n++)
			{
				BSPLable* label = dest->dest->labelSet[n];

				if (label->objvalue2 < rightbound && !hitlowerbound)
				{
					floatType lowerbound = __max(label->objvalue2, dDist->lambda_min);
					floatType tmpdmd = (dDist->GetCDF(rightbound) - dDist->GetCDF(lowerbound)) * dest->assDemand;
					ttmincost += label->objvalue1 * tmpdmd;
					rightbound = lowerbound;
					if (rightbound == dDist->lambda_min) hitlowerbound = true;
				}

			}
		}
	}

	RGapIndicator = fabs(1.0 - ttmincost / ttcost);
}

void	ENET::ConComputeRGII() // Outer loop
{
	EORIGIN* org;
	EDest* dest;
	DistanceDist* dDist;
	EPATH* path;
	OFV = 0; 
	numOfPath = 0;
	double numOfSubPath = 0.0;
	double numOfOD = 0.0;
	floatType   ttcost = 0.0;
	floatType	ttmincost = 0.0;
	for (int i = 0; i < numOfLink; i++)
	{
		OFV += linkVector[i]->GetIntCost();
	}
	double maxgap = -1;
	for (int i = 0; i < numOfOrg; i++)
	{ 
		org = orgVector[i];
		LabelSettingBSPTree(org);
		for (int j = 0; j < org->numOfDest; j++)
		{
			numOfOD += 1;
			dest = org->destVector[j];
			dest->UpdatePathSetCost();
			//if (dest->dest->id == 18 && org->org->id == 20&& curIter > 10)//
			//{
			//	cout << endl;
			//}
			floatType t1 = 0, t2 = 0, t3 = 0, t4 = 0;
			for (vector<EPATH*>::iterator it = dest->pathSet.begin(); it != dest->pathSet.end(); it++)
			{
				ttcost += (*it)->cost * (*it)->flow;
				t1 += (*it)->cost * (*it)->flow;
				t3 += (*it)->flow;
				//for (int i = 0; i < size((*it)->vLinks); i++)
				//{
					//if ((*it)->vLinks[i]->id == 880)
					//{
						//cout << "path id:"<< (*it)->id << "path distance:" << (*it)->length << "flow:" << (*it)->flow << endl;
					//}
				//}
			}
			// merge range
			for (int k = dest->rangeVector.size() - 1; k > 0; k--)
			{
				SubRange* last_range = dest->rangeVector[k];
				SubRange* lasttwo_range = dest->rangeVector[k-1];
				last_range->addpath = false;
				lasttwo_range->addpath = false;
				if (key_compare(last_range->GroupPaths,lasttwo_range->GroupPaths))
				{
					// set merge
					lasttwo_range->lowerBound = last_range->lowerBound;
					lasttwo_range->dmd += last_range->dmd;
					for (GroupPathMapIter it = lasttwo_range->GroupPaths.begin(); it != lasttwo_range->GroupPaths.end(); it++)
					{
						int pid = it->first;
						it->second->subpathflow += last_range->GroupPaths[pid]->subpathflow;
					}
					dest->rangeVector.erase(dest->rangeVector.begin() + k);
					//lasttwo_range->UpdateRangePathSet();
				}	
				else
				{
					//last_range->UpdateRangePathSet();
					
				}

			}
			dest->rangeVector[0]->addpath = false;
			
		/*	for (int k = dest->rangeVector.size() - 1; k >= 0; k--)
			{
				SubRange* range = dest->rangeVector[k];
				range->addpath = false;
			}*/
			
			//compute the minimum total cost 
			dDist = dest->dDist;
			floatType rightbound = dDist->lambda_max;
			bool hitlowerbound = false;
			sort(dest->dest->labelSet.begin(), dest->dest->labelSet.end(), lexicogr_order);
			
			for (int n = 0; n < dest->dest->labelSet.size(); n++)
			{
				BSPLable* label = dest->dest->labelSet[n];
				if (label->objvalue2 < rightbound && !hitlowerbound)
				{
					floatType lowerbound = __max(label->objvalue2, dDist->lambda_min);
					floatType tmpdmd = (dDist->GetCDF(rightbound) - dDist->GetCDF(lowerbound)) * dest->assDemand;
					ttmincost += label->objvalue1 * tmpdmd;
					t2 += label->objvalue1 * tmpdmd;
					rightbound = lowerbound;
					if (rightbound == dDist->lambda_min) hitlowerbound = true;
				}
				//if (dest->rangeVector.size() > 1)
				//{
				//	cout << endl;
				//}
				//column generation
				if (label->objvalue2 <= dDist->lambda_max)
				{
					for (int k = dest->rangeVector.size() - 1; k >= 0; k--)
					{	
						SubRange* range = dest->rangeVector[k];
						//add path to this range if it has not been included
						if (round(label->objvalue2*10)/10 <= round(10*range->lowerBound)/10 && !range->addpath)//&& label->objvalue1< range->GroupPaths[range->minIx]->path->cost
						{
							EPATH* newpath = InitializePathByLabel(dest->origin, dest, label);
							path = dest->InPathSet(newpath);

							if (path == NULL)
							{
								dest->pathSet.push_back(newpath);
								SubPath* newsubpath = new SubPath(newpath, 0.0);
								pair<GroupPathMapIter, bool> ret = range->GroupPaths.insert(pair<int, SubPath*>(newpath->id, newsubpath));
							}
							else
							{
								if (!range->Findpathbyid(path->id))
								{
									SubPath* subpath = new SubPath(path, 0.0);
									pair<GroupPathMapIter, bool> ret = range->GroupPaths.insert(pair<int, SubPath*>(path->id, subpath));
								}
							}
							range->addpath = true; // add path only once, which should be the lowest cost path
						}
						//set split
						if (round(label->objvalue2*10)/10 < round(range->upperBound*10)/10  && round(label->objvalue2*10)/10 > round(range->lowerBound*10)/10)
						{
							// geneate a new subrange and add this path to it
							floatType ndmd = (dDist->GetCDF(range->upperBound) - dDist->GetCDF(label->objvalue2)) * dest->assDemand;
							SubRange* newrange = new SubRange(label->objvalue2, range->upperBound, ndmd, dest);
							EPATH* newpath = InitializePathByLabel(dest->origin, dest, label);
							EPATH* tmpath = dest->InPathSet(newpath);
							if (!tmpath)// a new path
							{
								tmpath = newpath;
								dest->pathSet.push_back(tmpath);
							}
							SubPath* subpath = new SubPath(tmpath, 0.0);
							pair<GroupPathMapIter, bool> ret = newrange->GroupPaths.insert(pair<int, SubPath*>(tmpath->id, subpath));
							
							for (GroupPathMapIter gp = range->GroupPaths.begin(); gp != range->GroupPaths.end(); gp++)
							{
								SubPath* newsubpath = new SubPath(gp->second->path, ndmd / range->dmd * gp->second->subpathflow);
								pair<GroupPathMapIter, bool> ret = newrange->GroupPaths.insert(pair<int, SubPath*>(gp->second->path->id, newsubpath));
							}
							/*else
							{
								cout << "something wrong" << endl;
							}*/
							dest->rangeVector.insert(dest->rangeVector.begin() + k, newrange);
							newrange->addpath = true;
							range->upperBound = label->objvalue2;
							floatType dmd = (dDist->GetCDF(range->upperBound) - dDist->GetCDF(range->lowerBound)) * dest->assDemand;
							for (GroupPathMapIter gp = range->GroupPaths.begin(); gp != range->GroupPaths.end(); gp++)
							{
								gp->second->subpathflow = gp->second->subpathflow * dmd / range->dmd;
							}
							range->dmd = dmd;
						}

					}

				}
				if (label->objvalue2 <= dDist->lambda_min) break;
			
			}
			numOfSubPath += dest->rangeVector.size();
		}
	}
	double subRate = numOfSubPath / numOfOD;
	cout << "subNum:" << numOfSubPath << endl;
	cout << "odNum:" << numOfOD << endl;
	cout << "subNumRate:" << subRate << endl;
	RGapIndicator = fabs(1.0 - ttmincost / ttcost);
	//cout << "max gap of an OD:" << maxgap << endl;
}





void	ENET::ConComputeRG()
{
	EORIGIN* org;
	EDest* dest;
	DistanceDist* dDist;
	OFV = 0;
	numOfPath = 0;
	floatType   ttcost = 0.0;
	floatType	ttmincost = 0.0;
	EPATH* path;

	for (int i = 0; i < numOfLink; i++)
	{
		OFV += linkVector[i]->GetIntCost();
	}

	for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		LabelSettingBSPTree(org);
		for (int j = 0; j < org->numOfDest; j++)
		{
			dest = org->destVector[j];
			dest->UpdatePathSetCost();
			numOfPath += dest->pathSet.size();

			floatType t1 = 0, t2 = 0;
			for (vector<EPATH*>::iterator it = dest->pathSet.begin(); it != dest->pathSet.end(); it++)
			{
				ttcost += (*it)->cost * (*it)->flow;
				t1 += (*it)->cost * (*it)->flow;
				//cout << "path id:" << (*it)->id <<",len:"<< (*it)->length << ",flow:" << (*it)->flow << ",cost:"<< (*it)->cost<<endl;
			}

			dDist = dest->dDist;
			floatType rightbound = dDist->lambda_max;
			bool hitlowerbound = false;
			sort(dest->dest->labelSet.begin(), dest->dest->labelSet.end(), lexicogr_order);

			//if ( dest->dest->id == 100 && org->org->id == 3)
			//{
			////	dest->RangeMerge();
			//	cout << endl;
			//}

			for (int n = 0; n < dest->dest->labelSet.size(); n++)
			{
				BSPLable* label = dest->dest->labelSet[n];

				if (label->objvalue2 <= rightbound &&!hitlowerbound)
				{
					floatType lowerbound = __max(label->objvalue2, dDist->lambda_min);
					floatType tmpdmd = (dDist->GetCDF(rightbound) - dDist->GetCDF(lowerbound)) * dest->assDemand;
					ttmincost += label->objvalue1 * tmpdmd;
					t2 += label->objvalue1 * tmpdmd;

					rightbound = lowerbound;
					if (rightbound == dDist->lambda_min) hitlowerbound = true;
				}


			
				if (label->objvalue2 <= dDist->lambda_max && label->objvalue2 > dDist->lambda_min)
				{
					for(int k = dest->rangeVector.size()-1;k>=0;k--)
					{			
						SubRange* range = dest->rangeVector[k];
			
						//add path to this range if it has not been included
						if (label->objvalue2 <= range->lowerBound)
						{
							EPATH* newpath = InitializePathByLabel(dest->origin, dest, label);
							path = dest->InPathSet(newpath);

							if (path == NULL)
							{
								dest->pathSet.push_back(newpath);
								SubPath* newsubpath = new SubPath(newpath, 0.0);
								pair<GroupPathMapIter, bool> ret = range->GroupPaths.insert(pair<int, SubPath*>(newpath->id, newsubpath));
							}
							else
							{
								if (!range->Findpathbyid(path->id))
								{
									SubPath* subpath = new SubPath(path, 0.0);
									pair<GroupPathMapIter, bool> ret = range->GroupPaths.insert(pair<int, SubPath*>(path->id, subpath));
								}
							}				
						}

						//set split
						if (label->objvalue2 - range->upperBound < 0 && label->objvalue2 - range->lowerBound>0)
						{
							// geneate a new subrange and add this path to it
							floatType ndmd = (dDist->GetCDF(range->upperBound) - dDist->GetCDF(label->objvalue2)) * dest->assDemand;
							SubRange* newrange = new SubRange(label->objvalue2, range->upperBound, ndmd, dest);


							EPATH* newpath = InitializePathByLabel(dest->origin, dest, label);
							
							if (dest->InPathSet(newpath) == NULL)
							{

								SubPath* subpath = new SubPath(newpath, 0.0);
								pair<GroupPathMapIter, bool> ret = newrange->GroupPaths.insert(pair<int, SubPath*>(newpath->id, subpath));
								dest->pathSet.push_back(newpath);
								// add original subpaths to this new range
								if (newrange->dest->origin->org->id == 16 && newrange->dest->dest->id == 48)
								{
									cout << endl;
								}
								floatType sumflow = 0;
								for (GroupPathMapIter gp = range->GroupPaths.begin(); gp != range->GroupPaths.end(); gp++)
								{
									SubPath* newsubpath = new SubPath(gp->second->path, ndmd/range->dmd * gp->second->subpathflow);
									sumflow += newsubpath->subpathflow;
									pair<GroupPathMapIter, bool> ret = newrange->GroupPaths.insert(pair<int, SubPath*>(gp->second->path->id, newsubpath));
								}
								if (fabs(sumflow - ndmd) > 1e-8)
								{
									cout << "not conservative:" << fabs(sumflow - ndmd) << endl;
									cout << endl;
								}
							}
							else
							{
								//cout << "something wrong!" << endl;
							}
							dest->rangeVector.insert(dest->rangeVector.begin()+k,newrange);


							range->upperBound = label->objvalue2;
							floatType dmd = (dDist->GetCDF(range->upperBound) - dDist->GetCDF(range->lowerBound)) * dest->assDemand;
							for (GroupPathMapIter gp = range->GroupPaths.begin(); gp != range->GroupPaths.end(); gp++)
							{
								gp->second->subpathflow = gp->second->subpathflow * dmd / range->dmd;
							}
							range->dmd = dmd;
						}						
					}
				}

				if (label->objvalue2 <= dDist->lambda_min)
				{
					//add this new path to the last subrange if it is new
					for (int k = dest->rangeVector.size() - 1; k >= 0; k--)
					{
						SubRange* range = dest->rangeVector[k];
						EPATH* newpath = InitializePathByLabel(dest->origin, dest, label);
						path = dest->InPathSet(newpath);

						if (path == NULL)
						{
							dest->pathSet.push_back(newpath);
							SubPath* newsubpath = new SubPath(newpath, 0.0);
							pair<GroupPathMapIter, bool> ret = range->GroupPaths.insert(pair<int, SubPath*>(newpath->id, newsubpath));
						}
						else
						{
							if (!range->Findpathbyid(path->id))
							{
								SubPath* subpath = new SubPath(path, 0.0);
								pair<GroupPathMapIter, bool> ret = range->GroupPaths.insert(pair<int, SubPath*>(path->id, subpath));
							}
						}
					
					}
				}
				if (label->objvalue2 < dDist->lambda_min)
					break;


			}
			
			//dest->RangeMerge();

			//cout << "org:" << org->org->id << ",dest:" << dest->dest->id << ",t1:" << t1 << ",t2:" << t2 << endl << endl;

		}
	}

	RGapIndicator = fabs(1.0 - ttmincost / ttcost);

}

void EDest::RangeMerge()
{
	
	for (int k = rangeVector.size() - 1; k > 0; k--)
	{
		SubRange* range = rangeVector[k];
		SubRange* nextrange = rangeVector[k-1];
		if (range->GroupPaths == nextrange->GroupPaths)
		{
			cout << endl;
		}
		cout << endl;
	}
}


void	ENET::LabelCorrectingSP(EORIGIN* org)
{
	deque<ENODE*> nl;
	ENODE *node, *curNode, *scanNode;
	vector<ELINK*>::iterator pv;
	ENODE* rootNode = org->org;
	for (int j = 0; j < numOfNode; j++)
	{
		node = nodeVector[j];
		node->PathElem->cost = POS_INF_FLOAT;
		node->PathElem->via = NULL;
		node->scanStatus = 0;
	}

	rootNode->PathElem->cost = 0;
	nl.push_back(rootNode);

	while (!nl.empty())
	{
		curNode = nl.front();
		nl.pop_front(); //delete it from the deque;
		curNode->scanStatus = -1; //mark it as been used but not in queue right now.
		for (pv = curNode->forwStar.begin(); pv != curNode->forwStar.end(); pv++)
		{
			scanNode = (*pv)->head;
			floatType dp = curNode->PathElem->cost;
			floatType curTT = (*pv)->length;
			if (scanNode->PathElem->cost > curTT + dp)
			{
				scanNode->PathElem->cost = curTT + dp;
				scanNode->PathElem->via = *pv;

				if (scanNode->scanStatus == 0) // if never been used, insert it to the back of dq
				{
					nl.push_back(scanNode);
					scanNode->scanStatus = 1; // now being used
				}
				else if (scanNode->scanStatus == -1)
				{
					nl.push_front(scanNode);
					scanNode->scanStatus = 1; // now being used
				}
			}
		}
	}

}

void	ENET::IniRandomDistByFacor(floatType low_ratio, floatType up_ratio, DistanceDist::Dists type)
{

	EORIGIN* org;
	EDest* dest;
	for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		LabelCorrectingSP(org);
		for (int j = 0; j < org->numOfDest; j++)
		{
			dest = org->destVector[j];
			floatType min = dest->dest->PathElem->cost * low_ratio;
			floatType max = dest->dest->PathElem->cost * up_ratio;
			DistanceDist* dist = new DistanceDist(min, max, type);
			dest->dDist = dist;
		}
	}

}


void	ENET::IniRandomDistByLowUpValue(floatType low, floatType up, DistanceDist::Dists type)
{
	EORIGIN* org;
	EDest* dest;
	for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		for (int j = 0; j < org->numOfDest; j++)
		{
			dest = org->destVector[j];
			if (i == 0)
			{
				low = 25;
				up = 40;
			}
			else
			{

				low = 30;
				up = 60;
			}

			DistanceDist* dist = new DistanceDist(low, up, type);
			dest->dDist = dist;
		}
	}
}

void	ENET::InfinityDistance(floatType low, floatType up, DistanceDist::Dists type)
{
	EORIGIN* org;
	EDest* dest;
	for (int i = 0; i < numOfOrg; i++)
	{
		org = orgVector[i];
		for (int j = 0; j < org->numOfDest; j++)
		{
			dest = org->destVector[j];
			low = 1000000; 
			up = 10000000;
			DistanceDist* dist = new DistanceDist(low, up, type);
			dest->dDist = dist;
		}
	}
}


void	ENET::ReportIter()
{
	string IterConvName = networkName + "-" + GetAlgorithmName() + to_string(1.1) + ".conv";

	ofstream outfile;
	if (!TNM_OpenOutFile(outfile, IterConvName))
	{
		cout << "\n\tFail to prepare report: Cannot open .conv file to write convergence information!" << endl;
		return;
	}
	cout << "\tWriting iterative information into " << IterConvName << " for covergence info!" << endl;
	outfile << "Iter,time,RelativeGap,ofv,numOfPath,numOfSubpath,numOfTrips,innerIters" << endl;
	outfile << "0,0,1,0,0,0,0" << endl;
	//outfile << "0,1,,0,,,," << endl;
	for (int i = 0; i < itersRecord.size(); i++)
	{
		EITERELEM* it = itersRecord[i];
		outfile << TNM_IntFormat(it->iter) << "," << TNM_FloatFormat(it->time, 20, 18) << "," << TNM_FloatFormat(it->convRGap, 20, 18) << "," << TNM_FloatFormat(it->ofv, 20, 18)
			<< "," << TNM_IntFormat(it->numberofpaths) << "," << TNM_IntFormat(it->numberofsubpaths)
			<< "," << TNM_FloatFormat(it->numoftrips, 20, 18) << "," << TNM_IntFormat(it->innerIters) << endl;
	}
	outfile.close();
}

void	ENET::RecordCurrentIter()
{
	EITERELEM* iterElem = new EITERELEM;
	iterElem->iter = curIter;
	iterElem->time = 1.0 * (clock() - m_startRunTime) / CLOCKS_PER_SEC;
	iterElem->convRGap = RGapIndicator;
	iterElem->numberofpaths = numOfPath;
	iterElem->numberofsubpaths = numOfSubPath;
	iterElem->ofv = OFV;
	iterElem->numoftrips = numOfTrips;
	iterElem->innerIters = innerIters;
	itersRecord.push_back(iterElem);
}
