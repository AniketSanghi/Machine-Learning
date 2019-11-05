#include <iostream>
#include <fstream>
#include <string>

#include "timer.h"
#include "fastXML.h"

using namespace std;

void help()
{
	cerr<<"Sample Usage :"<<endl;
	cerr<<"./fastXML_train [feature file name] [label file name] [model folder name] -T 1 -s 0 -t 50 -b 1.0 -c 1.0 -m 10 -l 100 -q 0"<<endl<<endl;

	cerr<<"-T Number of threads to use. default=1"<<endl;
	cerr<<"-s Starting tree index. default=0"<<endl;
	cerr<<"-t Number of trees to be grown. default=50"<<endl;
	cerr<<"-b Feature bias value, extre feature value to be appended. default=1.0"<<endl;
	cerr<<"-c SVM weight co-efficient. default=1.0"<<endl;
	cerr<<"-m Maximum allowed instances in a leaf node. Larger nodes are attempted to be split, and on failure converted to leaves		default=10"<<endl;
	cerr<<"-l Number of label-probability pairs to retain in a leaf. default=100"<<endl;
	cerr<<"-q quiet option (0/1). default=0"<<endl;

	cerr<<"feature and label files are in sparse matrix format"<<endl;
	exit(1);
}

Param parse_param(_int argc, char* argv[])
{
	Param param;

	string opt;
	string sval;
	_float val;

	for(_int i=0; i<argc; i+=2)
	{
		opt = string(argv[i]);
		sval = string(argv[i+1]);
		val = stof(sval);

		if(opt=="-m")
			param.max_leaf = (_int)val;
		else if(opt=="-l")
			param.lbl_per_leaf = (_int)val;
		else if(opt=="-b")
			param.bias = (_float)val;
		else if(opt=="-c")
			param.log_loss_coeff = (_float)val;
		else if(opt=="-T")
			param.num_thread = (_int)val;
		else if(opt=="-s")
			param.start_tree = (_int)val;
		else if(opt=="-t")
			param.num_tree = (_int)val;
		else if(opt=="-q")
			param.quiet = (_bool)val;
	}

	return param;
}

int main(int argc, char* argv[])
{
	if(argc < 4)
		help();

	string ft_file = string(argv[1]);
	check_valid_filename(ft_file, true);
	SMatF* trn_X_Xf = new SMatF(ft_file);

	string lbl_file = string(argv[2]);
	check_valid_filename(lbl_file, true);
	SMatF* trn_X_Y = new SMatF(lbl_file);

	string model_folder = string(argv[3]);
	check_valid_foldername(model_folder);

	string prop_file = string(argv[4]);
	check_valid_filename(prop_file, true);
	ifstream fin;
	fin.open(prop_file);
	VecF inv_props;
	for(_int i=0; i<trn_X_Y->nr; i++)
	{
		_float f;
		fin>>f;
		inv_props.push_back(f);
	}
	fin.close();

	Param param = parse_param(argc-5,argv+5);
	param.num_Xf = trn_X_Xf->nr;
	param.num_Y = trn_X_Y->nr;
	param.write(model_folder+"/param");

	


	for(_int i=0; i<trn_X_Y->nc; i++)
		for(_int j=0; j<trn_X_Y->size[i]; j++)
			trn_X_Y->data[i][j].second *= inv_props[trn_X_Y->data[i][j].first];





	if( param.quiet )
		loglvl = LOGLVL::QUIET;

	_float train_time;
	train_trees( trn_X_Xf, trn_X_Y, param, model_folder, train_time );
	cout << "training time: " << train_time/3600.0 << " hr" << endl;

	delete trn_X_Xf;
	delete trn_X_Y;
}
