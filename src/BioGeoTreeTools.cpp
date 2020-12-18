/*
 * PhyloTree.cpp
 *
 *  Created on: Aug 15, 2009
 *      Author: Stephen A. Smith
 */

#include "BioGeoTreeTools.h"
#include "RateMatrixUtils.h"
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cmath>
using namespace std;

#include "superdouble.h"
#include "tree.h"
#include "tree_reader.h"
#include "node.h"
#include "vector_node_object.h"
#include "string_node_object.h"
#include "Utils.h"

Tree * BioGeoTreeTools::getTreeFromString(string treestring){
    TreeReader tr;
    return tr.readTree(treestring);
}

vector<Node *> BioGeoTreeTools::getAncestors(Tree & tree, Node & nodeId){
    vector<Node *> nodes;
    Node * current = &nodeId;
    while(current->hasParent()){
	current = current->getParent();
	nodes.push_back(current);
    }
    return nodes;
}

void BioGeoTreeTools::summarizeSplits(Node * node,map<vector<int>,vector<AncSplit> > & ans,map<int,string> &areanamemaprev, RateModel * rm){
    Superdouble best(0);
    Superdouble sum(0);
    map<Superdouble,string > printstring;
    int areasize = (*ans.begin()).first.size();
    map<int, vector<int> > * distmap = rm->get_int_dists_map(); 
    vector<int> bestldist;
    vector<int> bestrdist;
    map<vector<int>,vector<AncSplit> >::iterator it;
    bool first = true;
    for(it=ans.begin();it!=ans.end();it++){
	vector<int> dis = (*it).first;
	vector<AncSplit> tans = (*it).second;
	for (unsigned int i=0;i<tans.size();i++){
	    if (first == true){
		first = false;
		best = tans[i].getLikelihood();
		bestldist = (*distmap)[tans[i].ldescdistint];//tans[i].getLDescDist();
		bestrdist = (*distmap)[tans[i].rdescdistint];//tans[i].getRDescDist();
	    }else if (tans[i].getLikelihood() > best){
		best = tans[i].getLikelihood();
		bestldist = (*distmap)[tans[i].ldescdistint];//tans[i].getLDescDist();
		bestrdist = (*distmap)[tans[i].rdescdistint];//tans[i].getRDescDist();
	    }
	    //cout << -log(tans[i].getLikelihood()) << endl;
	    sum += tans[i].getLikelihood();
	}
    }
    Superdouble test2(2);
    for(it=ans.begin();it!=ans.end();it++){
	vector<int> dis = (*it).first;
	vector<AncSplit> tans = (*it).second;
	for (unsigned int i=0;i<tans.size();i++){
	    if ((best.getLn()-(tans[i].getLikelihood().getLn())) < test2){
		string tdisstring ="";
		int  count1 = 0;
		for(int m=0;m<areasize;
		    //tans[i].getLDescDist().size();
		    m++){
		    //if(tans[i].getLDescDist()[m] == 1){
		    if((*distmap)[tans[i].ldescdistint][m] == 1){
			tdisstring += areanamemaprev[m];
			count1 += 1;
			if(count1 < calculate_vector_int_sum(&(*distmap)[tans[i].ldescdistint])){
			    tdisstring += "_";
			}
		    }

		}tdisstring += "|";
		count1 = 0;
		for(int m=0;m<areasize;
		    //tans[i].getRDescDist().size();
		    m++){
		    if((*distmap)[tans[i].rdescdistint][m] == 1){
			tdisstring += areanamemaprev[m];
			count1 += 1;
			if(count1 < calculate_vector_int_sum(&(*distmap)[tans[i].rdescdistint])){
			    tdisstring += "_";
			}
		    }
		}
		printstring[tans[i].getLikelihood()] = tdisstring;
	    }
	}
    }
    Superdouble none(-1);
    map<Superdouble,string >::reverse_iterator pit;
    for(pit=printstring.rbegin();pit != printstring.rend();pit++){
	Superdouble lnl(((*pit).first));
	cout << "\t" << (*pit).second << "\t" << double(lnl/sum) << "\t(" << double(none*lnl.getLn())<< ")"<< endl;
    }
    StringNodeObject disstring ="";
    int  count = 0;
    for(unsigned int m=0;m<bestldist.size();m++){
	if(bestldist[m] == 1){
	    disstring += areanamemaprev[m];
	    count += 1;
	    //if(count < calculate_vector_int_sum(&bestldist))
	    if(count < accumulate(bestldist.begin(),bestldist.end(),0))
		disstring += "_";
	}
    }disstring += "|";
    count = 0;
    for(unsigned int m=0;m<bestrdist.size();m++){
	if(bestrdist[m] == 1){
	    disstring += areanamemaprev[m];
	    count += 1;
	    //if(count < calculate_vector_int_sum(&bestrdist))
	    if(count < accumulate(bestrdist.begin(),bestrdist.end(),0))
		disstring += "_";
	}
    }
    string spl = "split";
    node->assocObject(spl,disstring);
    //cout << -log(best) << " "<< best/sum << endl;
}


void BioGeoTreeTools::summarizeAncState(Node * node,vector<Superdouble> & ans,map<int,string> &areanamemaprev, RateModel * rm, nlohmann::json& output){
    //Superdouble best(ans[1]);
    Superdouble best(ans[1]);//use ans[1] because ans[0] is just 0
    Superdouble sum(0);
    map<Superdouble,string > printstring;
    int areasize = rm->get_num_areas();
    map<int, vector<int> > * distmap = rm->get_int_dists_map(); 
    vector<int> bestancdist;
    Superdouble zero(0);
    for(unsigned int i=1;i<ans.size();i++){//1 because 0 is just the 0 all extinct one
	if (ans[i] >= best && ans[i] != zero){ //added != 0, need to test
	    best = ans[i];
	    bestancdist = (*distmap)[i];
	}
	sum += ans[i];
    }
    for(unsigned int i=1;i<ans.size();i++){
	nlohmann::json tmp_json;
	tmp_json["llh"] = double(ans[i].getLn());
	tmp_json["ratio"] = double(ans[i]/sum);
	tmp_json["distribution"] = dist_to_int(distmap->at(i));
	output.push_back(tmp_json);
    }
    Superdouble none(-1);
    Superdouble test2(2);
    for(unsigned int i=0;i<ans.size();i++){
	if (((best.getLn())-(ans[i].getLn()) ) < test2){
	    string tdisstring ="";
	    int  count1 = 0;
	    for(int m=0;m<areasize;m++){
		if((*distmap)[i][m] == 1){
		    tdisstring += areanamemaprev[m];
		    count1 += 1;
		    if(count1 < calculate_vector_int_sum(&(*distmap)[i])){
			tdisstring += "_";
		    }
		}
	    }
	    printstring[ans[i]] = tdisstring;
	}
    }

    map<Superdouble,string >::reverse_iterator pit;
    for(pit=printstring.rbegin();pit != printstring.rend();pit++){
	Superdouble lnl(((*pit).first));
	//cout << lnl << endl;
	cout << "\t" << (*pit).second << "\t" << double(lnl/sum) << "\t(" << double(none*lnl.getLn()) << ")"<< endl;
    }
    StringNodeObject disstring ="";
    int  count = 0;
    for(unsigned int m=0;m<bestancdist.size();m++){
	if(bestancdist[m] == 1){
	    disstring += areanamemaprev[m];
	    count += 1;
	    //if(count < calculate_vector_int_sum(&bestldist))
	    if(count < accumulate(bestancdist.begin(),bestancdist.end(),0))
		disstring += "_";
	}
    }
    string spl = "state";
    node->assocObject(spl,disstring);
    //cout << -log(best) << " "<< best/sum << endl;	
}

string BioGeoTreeTools::get_string_from_dist_int(int dist,map<int,string> &areanamemaprev, RateModel * rm){
    map<int, vector<int> > * distmap = rm->get_int_dists_map();
    vector<int> bestancdist = (*distmap)[dist];

    StringNodeObject disstring ="";
    int  count = 0;
    for(unsigned int m=0;m<bestancdist.size();m++){
	if(bestancdist[m] == 1){
	    disstring += areanamemaprev[m];
	    count += 1;
	    //if(count < calculate_vector_int_sum(&bestldist))
	    if(count < accumulate(bestancdist.begin(),bestancdist.end(),0))
		disstring += "_";
	}
    }
    return disstring;
}

