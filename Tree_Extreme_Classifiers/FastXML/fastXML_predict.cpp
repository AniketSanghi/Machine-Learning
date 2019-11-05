#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <map>
#include <random>
#include <thread>
#include <mutex>
#include <functional>
#include <cassert>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <ctime>

#include <chrono>

// Config.h
	#define _int int
	#define _float float
	#define _llint long long int
	#define _double double
	#define _bool bool


	using namespace std;
	using namespace std::chrono;


// util.h
		typedef pair<_int,_float> pairIF;
		typedef pair<_int,_double> pairID;
		typedef pair<_int,_int> pairII;
		typedef pair<_int,_bool> pairIB;

		#define LINE (cout<<__LINE__<<endl)
		#define NAME_LEN 1000
		#define SQ(x) ((x)*(x))
		#define INF FLT_MAX
		#define NEG_INF FLT_MIN
		#define EPS 1e-10

		enum LOGLVL {QUIET, PROGRESS, DEBUG};

		template <typename T1,typename T2>
		_bool comp_pair_by_second_desc(pair<T1,T2> a, pair<T1,T2> b)
		{
			if(a.second>b.second)
				return true;
			return false;
		}

		template <typename T1,typename T2>
		_bool comp_pair_by_second(pair<T1,T2> a, pair<T1,T2> b)
		{
			if(a.second<b.second)
				return true;
			return false;
		}

		template <typename T1,typename T2>
		_bool comp_pair_by_first(pair<T1,T2> a, pair<T1,T2> b)
		{
			if(a.first<b.first)
				return true;
			return false;
		}

		template <typename T>
		void Realloc(_int old_size, _int new_size, T*& vec)
		{
			T* new_vec = new T[new_size];
			_int size = min(old_size,new_size);
			copy_n(vec,size,new_vec);
			delete [] vec;
			vec = new_vec;
		}

		template <typename T>
		void copy_S_to_D(_int size, pair<_int,T>* sarr, T* darr)
		{
			for(_int i=0; i<size; i++)
			{
				darr[sarr[i].first] = sarr[i].second;
			}
		}

		template <typename T>
		void reset_D(_int size, pair<_int,T>* sarr, T* darr)
		{
			for(_int i=0; i<size; i++)
			{
				darr[sarr[i].first] = 0;
			}
		}

		inline void check_valid_filename(string fname, _bool read=true)
		{
			_bool valid;
			ifstream fin;
			ofstream fout;
			if(read)
			{
				fin.open(fname);
				valid = fin.good();
			}
			else
			{
				fout.open(fname);
				valid = fout.good();
			}

			if(!valid)
			{
				cerr<<"error: invalid file name: "<<fname<<endl<<"exiting..."<<endl;
				exit(1);
			}
			if(read)
			{
				fin.close();
			}
			else
			{
				fout.close();
			}
		}

		inline void check_valid_foldername(string fname)
		{
			string tmp_file = fname+"/tmp.txt";
			ofstream fout(tmp_file);

			if(!fout.good())
			{
				cerr<<"error: invalid folder name: "<<fname<<endl<<"exiting..."<<endl;
				exit(1);
			}
			remove(tmp_file.c_str());
		}

		template <typename T>
		void print_vector( ostream& fout, vector<T>& vec )
		{
			for( _int i=0; i<vec.size(); i++ )
				if(i==0)
					fout << vec[i];
				else
					fout << " " << vec[i];
			fout << endl;
		}

		template <typename T1, typename T2>
		void print_vector( ostream& fout, vector< pair<T1,T2> >& vec )
		{
			for( _int i=0; i<vec.size(); i++ )
				if(i==0)
					fout << vec[i].first << ":" << vec[i].second;
				else
					fout << " " << vec[i].first << ":" << vec[i].second;
			fout << endl;
		}

// timer.h

		class Timer
			{
				private:
					high_resolution_clock::time_point start_time;
					high_resolution_clock::time_point stop_time;

				public:
					Timer()
					{
					}

					void tic()
					{
						start_time = high_resolution_clock::now();
					}

					_float toc()
					{
						stop_time = high_resolution_clock::now();
						_float elapsed_time = duration_cast<microseconds>(stop_time - start_time).count() / 1000000.0;
						return elapsed_time;
					}
			};



// mat.h

			typedef vector<_int> VecI;
		typedef vector<_float> VecF;
		typedef vector<_double> VecD;
		typedef vector<pairII> VecII;
		typedef vector<pairIF> VecIF;
		typedef vector<_bool> VecB;

		template <typename T>
		class SMat // a column-major sparse matrix of type T
		{
		public:
			_bool contiguous;
			_int nc;
			_int nr;
			_int* size;
			pair<_int,T>** data;
			pair<_int,T>* cdata;

			SMat( _bool contiguous = false )
			{
				this->contiguous = contiguous;
				nc = 0;
				nr = 0;
				size = NULL;
				data = NULL;
				cdata = NULL;
			}

			SMat( _int nr, _int nc, _bool contiguous = false )
			{
				this->contiguous = contiguous;
				this->nr = nr;
				this->nc = nc;
				size = new _int[nc]();
				data = new pair<_int,T>*[nc];
				for(_int i=0; i<nc; i++)
					data[i] = NULL;
				cdata = NULL;
			}

			SMat( _int nr, _int nc, _int nnz, _bool contiguous = false )
			{
				this->contiguous = contiguous;
				this->nr = nr;
				this->nc = nc;
				size = new _int[nc]();
				data = new pair<_int,T>*[nc];
				for(_int i=0; i<nc; i++)
					data[i] = NULL;

				if( contiguous )
					cdata = new pair<_int,T>[ nnz ];
				else
					cdata = NULL;
			}
			
			SMat(SMat<T>* mat)
			{
				contiguous = false;
				nc = mat->nc;
				nr = mat->nr;
				size = new _int[nc];

				for(_int i=0; i<nc; i++)
					size[i] = mat->size[i];

				data = new pair<_int,T>*[nc];
				for(_int i=0; i<nc; i++)
				{
					data[i] = new pair<_int,T>[size[i]];
					for(_int j=0; j<size[i]; j++)
					{
						data[i][j] = mat->data[i][j];
					}	
				}
			}

			friend istream& operator>>( istream& fin, SMat<T>& mat )
			{
				vector<_int> inds;
				vector<T> vals;

				fin >> mat.nc >> mat.nr;

				mat.size = new _int[ mat.nc ];
				mat.data = new pair<_int,T>*[ mat.nc ];

				fin.ignore();
				for(_int i=0; i<mat.nc; i++)
				{
					inds.clear();
					vals.clear();
					string line;
					getline(fin,line);
					line += "\n";
					_int pos = 0;
					_int next_pos;

					while(next_pos=line.find_first_of(": \n",pos))
					{
						if((size_t)next_pos==string::npos)
							break;
						inds.push_back(stoi(line.substr(pos,next_pos-pos)));
						pos = next_pos+1;

						next_pos = line.find_first_of(": \n",pos);
						if((size_t)next_pos==string::npos)
							break;

						vals.push_back(stof(line.substr(pos,next_pos-pos)));
						pos = next_pos+1;

					}

					assert(inds.size()==vals.size());
					assert(inds.size()==0 || inds[inds.size()-1]<mat.nr);

					mat.size[i] = inds.size();
					mat.data[i] = new pair<_int,T>[inds.size()];

					for(_int j=0; j<mat.size[i]; j++)
					{
						mat.data[i][j].first = inds[j];
						mat.data[i][j].second = (T)vals[j];
					}
				}

				return fin;
			}

			SMat(string fname)
			{
				contiguous = false;
				check_valid_filename(fname,true);

				ifstream fin;
				fin.open(fname);		

				vector<_int> inds;
				vector<T> vals;

				string lineTemp;
				getline(fin,lineTemp);
				getline(fin,lineTemp);
				getline(fin,lineTemp);

				char c;
				fin>>c;

				fin>>nc;
				fin>>nr;
				size = new _int[nc];
				data = new pair<_int,T>*[nc];

				fin.ignore();
				for(_int i=0; i<nc; i++)
				{
					inds.clear();
					vals.clear();
					string line;
					getline(fin,line);
					line += "\n";
					line.erase(0,1);
					_int pos = 0;
					_int next_pos;

					while(next_pos=line.find_first_of(": \n",pos))
					{
						if((size_t)next_pos==string::npos)
							break;
						inds.push_back(stoi(line.substr(pos,next_pos-pos)));
						pos = next_pos+1;

						next_pos = line.find_first_of(": \n",pos);
						if((size_t)next_pos==string::npos)
							break;

						vals.push_back(stof(line.substr(pos,next_pos-pos)));
						pos = next_pos+1;

					}

					assert(inds.size()==vals.size());
					assert(inds.size()==0 || inds[inds.size()-1]<nr);

					size[i] = inds.size();
					data[i] = new pair<_int,T>[inds.size()];

					for(_int j=0; j<size[i]; j++)
					{
						data[i][j].first = inds[j];
						data[i][j].second = (T)vals[j];
					}
				}

				fin.close();
			}

			_float get_ram()
			{
				_float ram = 0;
				ram += 4*2 + 8*2 + 4*nc + 8*nc;

				for( _int i=0; i<nc; i++ )
					ram += 8*size[i];

				return ram;
			}

			SMat<T>* transpose()
			{
				SMat<T>* tmat = new SMat<T>;
				tmat->nr = nc;
				tmat->nc = nr;
				tmat->size = new _int[tmat->nc]();
				tmat->data = new pair<_int,T>*[tmat->nc];

				for(_int i=0; i<nc; i++)
				{
					for(_int j=0; j<size[i]; j++)
					{
						tmat->size[data[i][j].first]++;
					}
				}

				for(_int i=0; i<tmat->nc; i++)
				{
					tmat->data[i] = new pair<_int,T>[tmat->size[i]];
				}

				_int* count = new _int[tmat->nc]();
				for(_int i=0; i<nc; i++)
				{
					for(_int j=0; j<size[i]; j++)
					{
						_int ind = data[i][j].first;
						T val = data[i][j].second;

						tmat->data[ind][count[ind]].first = i;
						tmat->data[ind][count[ind]].second = val;
						count[ind]++;
					}
				}

				delete [] count;
				return tmat;
			}

			void threshold( _float th )
			{
				for( _int i=0; i<nc; i++ )
				{
					_int count = 0;
					for( _int j=0; j<size[i]; j++ )
						count += fabs( data[i][j].second )>th;

					pair<_int,T>* newvec = new pair<_int,T>[count];
					count = 0;
					for( _int j=0; j<size[i]; j++ )
					{
						_int id = data[i][j].first;
						T val = data[i][j].second;
						if( fabs(val)>th )
							newvec[ count++ ] = make_pair( id, val );
					}
					size[i] = count;
					delete [] data[i];
					data[i] = newvec;
				}
			}

			void unit_normalize_columns()
			{
				for(_int i=0; i<nc; i++)
				{
					T normsq = 0;
					for(_int j=0; j<size[i]; j++)
						normsq += SQ(data[i][j].second);
					normsq = sqrt(normsq);

					if(normsq==0)
						normsq = 1;

					for(_int j=0; j<size[i]; j++)
						data[i][j].second /= normsq;
				}
			}

			vector<T> column_norms()
			{
				vector<T> norms(nc,0);

				for(_int i=0; i<nc; i++)
				{
					T normsq = 0;
					for(_int j=0; j<size[i]; j++)
						normsq += SQ(data[i][j].second);
					norms[i] = sqrt(normsq);
				}

				return norms;
			}

			~SMat()
			{
				delete [] size;

				if( contiguous )
				{
					delete [] cdata;
					delete [] data; 
				}
				else
				{
					for( _int i=0; i<nc; i++ )
						delete [] data[i];
					delete [] data;
				}
			}

			friend ostream& operator<<( ostream& fout, const SMat<T>& mat )
			{
				_int nc = mat.nc;
				_int nr = mat.nr;
				_int* size = mat.size;
				pairIF** data = mat.data;

				fout<<nc<<" "<<nr<<endl;

				fout << fixed << setprecision( 6 );
				for(_int i=0; i<nc; i++)
				{
					for(_int j=0; j<size[i]; j++)
					{
						if(j==0)
							fout<<data[i][j].first<<":"<<data[i][j].second;
						else
							fout<<" "<<data[i][j].first<<":"<<data[i][j].second;
					}
					fout<<endl;
				}

				return fout;
			}

			void write(string fname)
			{
				check_valid_filename(fname,false);

				ofstream fout;
				fout.open(fname);

				fout << (*this);

				fout.close();
			}

			void add(SMat<T>* smat)
			{
				if(nc != smat->nc || nr != smat->nr)
				{
					cerr<<"SMat::add : Matrix dimensions do not match"<<endl;
					cerr<<"Matrix 1: "<<nc<<" x "<<nr<<endl;
					cerr<<"Matrix 2: "<<smat->nc<<" x "<<smat->nr<<endl;
					exit(1);
				}

				bool* ind_mask = new bool[nr]();
				T* sum = new T[nr]();

				for(_int i=0; i<nc; i++)
				{
					vector<_int> inds;
					for(_int j=0; j<size[i]; j++)
					{
						_int ind = data[i][j].first;
						T val = data[i][j].second;

						sum[ind] += val;
						if(!ind_mask[ind])
						{
							ind_mask[ind] = true;
							inds.push_back(ind);
						}
					}

					for(_int j=0; j<smat->size[i]; j++)
					{
						_int ind = smat->data[i][j].first;
						T val = smat->data[i][j].second;

						sum[ind] += val;
						if(!ind_mask[ind])
						{
							ind_mask[ind] = true;
							inds.push_back(ind);
						}
					}

					sort(inds.begin(), inds.end());
					Realloc(size[i], inds.size(), data[i]);
					for(_int j=0; j<inds.size(); j++)
					{
						_int ind = inds[j];
						data[i][j] = make_pair(ind,sum[ind]);
						ind_mask[ind] = false;
						sum[ind] = 0;
					}
					size[i] = inds.size();
				}

				delete [] ind_mask;
				delete [] sum;
			}

			SMat<T>* prod(SMat<T>* mat2)
			{
				_int dim1 = nr;
				_int dim2 = mat2->nc;

				assert(nc==mat2->nr);

				SMat<T>* prodmat = new SMat<T>(dim1,dim2);
				vector<T> sum(dim1,0);

				for(_int i=0; i<dim2; i++)
				{
					vector<_int> indices;
					for(_int j=0; j<mat2->size[i]; j++)
					{
						_int ind = mat2->data[i][j].first;
						T prodval = mat2->data[i][j].second;

						for(_int k=0; k<size[ind]; k++)
						{
							_int id = data[ind][k].first;
							T val = data[ind][k].second;

							if(sum[id]==0)
								indices.push_back(id);

							sum[id] += val*prodval;
						}
					}

					sort(indices.begin(), indices.end());

					_int siz = indices.size();
					prodmat->size[i] = siz;
					prodmat->data[i] = new pair<_int,T>[siz];

					for(_int j=0; j<indices.size(); j++)
					{
						_int id = indices[j];
						T val = sum[id];
						prodmat->data[i][j] = make_pair(id,val);
						sum[id] = 0;
					}
				}

				return prodmat;
			}

			void append_bias_feat( _float bias_feat ) // assumes non-contiguous matrix
			{
				for( _int i=0; i<nc; i++ )
				{
					_int siz = size[i];
					Realloc( siz, siz+1, data[i] );
					data[i][siz] = make_pair( nr, bias_feat );
					size[i]++;	
				}
				nr++;
			}

			void active_dims( VecI& cols, VecI& dims, VecI& counts, VecI& countmap )
			{
				dims.clear();
				counts.clear();

				for( _int i=0; i<cols.size(); i++ )
				{
					_int inst = cols[i];
					for( _int j=0; j<size[inst]; j++ )
					{
						_int dim = data[inst][j].first;
						if( countmap[ dim ]==0 )
							dims.push_back(dim);
						countmap[ dim ]++;
					}
				}

				sort(dims.begin(),dims.end());

				for( _int i=0; i<dims.size(); i++ )
				{
					counts.push_back( countmap[ dims[i] ] );
					countmap[ dims[i] ] = 0;
				}
			}
		    
		    void shrink_mat( VecI& cols, SMat<T>*& s_mat, VecI& rows, VecI& countmap, _bool transpose )
		    {
		        _int s_nc = cols.size();
		        VecI counts;
		        active_dims( cols, rows, counts, countmap );

		        _int nnz = 0;
		        for( _int i=0; i<counts.size(); i++ )
		            nnz += counts[i];

		        _int* maps = new _int[ nr ];
		        for( _int i=0; i<rows.size(); i++ )
		            maps[ rows[i] ] = i;

		        _int s_nr = rows.size();
		        
		        if( transpose )
		        {
		            s_mat = new SMat<T>( s_nc, s_nr, nnz, true );
		        
		            _int sumsize = 0;
		            for( _int i=0; i<s_nr; i++ )
		            {
		                s_mat->size[i] = counts[i];
		                s_mat->data[i] = s_mat->cdata + sumsize;
		                sumsize += counts[i];
		            }
		            
		            for( _int i=0; i<s_nr; i++ )
		                counts[i] = 0;
		        }
		        else
		        {
		            s_mat = new SMat<T>( s_nr, s_nc, nnz, true );

		            _int sumsize = 0;
		            for( _int i=0; i<s_nc; i++)
		            {
		                _int col = cols[i];
		                s_mat->size[i] = size[ col ];
		                s_mat->data[i] = s_mat->cdata + sumsize;
		                sumsize += size[ col ];
		            }
		        }
		            
		        for( _int i=0; i<s_nc; i++ )
		        {	
		            _int col = cols[ i ];
		            for( _int j=0; j<size[ col ]; j++ )
		            {
		                _int row = maps[ data[ col ][ j ].first ];
		                _float val = data[ col ][ j ].second;
		                
		                if( transpose )
		                {
		                    s_mat->data[row][counts[row]] = make_pair( i, val );
		                    counts[row]++;
		                }
		                else
		                    s_mat->data[i][j] = make_pair( row, val );
		            }
		        }

		        delete [] maps;
		    }
		};
		                  
		template <typename T>
		class DMat // a column-major dense matrix of type T
		{
		public:
			_int nc;
			_int nr;
			T** data;

			DMat()
			{
				nc = 0;
				nr = 0;
				data = NULL;
			}

			DMat(_int nc, _int nr)
			{
				this->nc = nc;
				this->nr = nr;
				data = new T*[nc];
				for(_int i=0; i<nc; i++)
					data[i] = new T[nr]();
			}

			DMat(SMat<T>* mat)
			{
				nc = mat->nc;
				nr = mat->nr;
				data = new T*[nc];
				for(_int i=0; i<nc; i++)
					data[i] = new T[nr]();

				for(_int i=0; i<mat->nc; i++)
				{
					pair<_int,T>* vec = mat->data[i];
					for(_int j=0; j<mat->size[i]; j++)
					{
						data[i][vec[j].first] = vec[j].second;
					}
				}
			}

			~DMat()
			{
				for(_int i=0; i<nc; i++)
					delete [] data[i];
				delete [] data;
			}
		};

		typedef SMat<_float> SMatF;
		typedef DMat<_float> DMatF;

		void reindex_VecIF( VecIF& vec, VecI& index );

		void reindex_VecIF( VecIF& vec, VecI& index )
			{
			    for( _int i=0; i<vec.size(); i++ )
			        vec[i].first = index[ vec[i].first ];
			    return;
			}


// fastXML.h

			extern LOGLVL loglvl;
			extern mutex mtx;
			extern _bool USE_IDCG;
			extern thread_local mt19937 reng; // random number generator used during training 
			extern thread_local VecF discounts;
			extern thread_local VecF csum_discounts;
			extern thread_local VecF dense_w;
			extern thread_local VecI countmap;

			class Param
			{
			public:
				_int num_Xf;
				_int num_Y;
				_float log_loss_coeff;
				_int max_leaf;
				_int lbl_per_leaf;
				_float bias;
				_int num_thread;
				_int start_tree;
				_int num_tree;
				_bool quiet;

				Param()
				{
					num_Xf = 0;
					num_Y = 0;
					log_loss_coeff = 1.0;
					max_leaf = 10;
					lbl_per_leaf = 100;
					bias = 1.0;
					num_thread = 1;
					start_tree = 0;
					num_tree = 50;
					quiet = false;
				}

				Param(string fname)
				{
					check_valid_filename(fname,true);
					ifstream fin;
					fin.open(fname);
					fin >> (*this);
					fin.close();
				}

				void write(string fname)
				{
					check_valid_filename(fname,false);
					ofstream fout;
					fout.open(fname);
					fout << (*this); 
					fout.close();
				}

				friend istream& operator>>( istream& fin, Param& param )
				{
					fin >> param.num_Xf;
					fin >> param.num_Y;
					fin >> param.log_loss_coeff;
					fin >> param.max_leaf;
					fin >> param.lbl_per_leaf;
					fin >> param.bias;
					fin >> param.num_thread;
					fin >> param.start_tree;
					fin >> param.num_tree;
					fin >> param.quiet;
					return fin;
				}

				friend ostream& operator<<( ostream& fout, const Param& param )
				{
					fout << param.num_Xf << "\n";
					fout << param.num_Y << "\n";
					fout << param.log_loss_coeff << "\n";
					fout << param.max_leaf << "\n";
					fout << param.lbl_per_leaf << "\n";
					fout << param.bias << "\n";
					fout << param.num_thread << "\n";
					fout << param.start_tree<< "\n";
					fout << param.num_tree << "\n";
					fout << param.quiet << endl;
					return fout;
				}
			};

			class Node
			{
			public:
				_bool is_leaf;
				_int pos_child;
				_int neg_child;
				_int depth;
				VecI X;
				VecIF w;
				VecIF leaf_dist;

				Node()
				{
					is_leaf = false;
					depth = 0;
					pos_child = neg_child = -1;
				}

				Node(VecI X, _int depth, _int max_leaf)
				{
					this->X = X;
					this->depth = depth;
					this->pos_child = -1;
					this->neg_child = -1;

					if(X.size()<=max_leaf)
						this->is_leaf = true;
					else
						this->is_leaf = false;
				}

				~Node()
				{
				}

				_float get_ram()
				{
					_float ram = sizeof( Node );
					ram += X.size() * sizeof( _int );
					ram += w.size() * sizeof( pairIF );
					ram += leaf_dist.size() * sizeof( pairIF );
					return ram;
				}

				friend ostream& operator<<(ostream& fout, const Node& node)
				{
					fout<<(node.is_leaf?1:0)<<"\n";

					fout<<node.pos_child<<" "<<node.neg_child<<"\n";
					fout<<node.depth<<"\n";

					fout<<node.X.size();
					for(_int i=0; i<node.X.size(); i++)
						fout<<" "<<node.X[i];
					fout<<"\n";

					if(node.is_leaf)
					{
						fout<<node.leaf_dist.size();
						for(_int i=0; i<node.leaf_dist.size(); i++)
						{
							fout<<" "<<node.leaf_dist[i].first<<":"<<node.leaf_dist[i].second;
						}
						fout<<"\n";
					}
					else
					{
						fout<<node.w.size();
						for(_int i=0; i<node.w.size(); i++)
						{
							fout<<" "<<node.w[i].first<<":"<<node.w[i].second;
						}
						fout << "\n";
					}
					return fout;
				}

				friend istream& operator>>(istream& fin, Node& node)
				{
					fin>>node.is_leaf;
					fin>>node.pos_child>>node.neg_child>>node.depth;

					_int siz;
					_int ind;
					_float val;
					char c;

					node.X.clear();
					fin>>siz;
					for(_int i=0; i<siz; i++)
					{
						fin>>ind;	
						node.X.push_back(ind);
					}

					if(node.is_leaf)
					{
						node.leaf_dist.clear();
						fin>>siz;
						for(_int i=0; i<siz; i++)
						{
							fin>>ind>>c>>val;
							node.leaf_dist.push_back(make_pair(ind,val));
						}
					}
					else
					{
						node.w.clear();
						fin>>siz;
						for(_int i=0; i<siz; i++)
						{
							fin>>ind>>c>>val;
							node.w.push_back(make_pair(ind,val));
						}	
					}
					return fin;
				}
			};

			template <typename GNode>
			class GTree // General instance tree with any kind of GNode. Supports FastXML, PfastreXML and SwiftXML
			{
			public:
				_int num_Xf;
				_int num_Y;
				vector<GNode*> nodes;

				GTree()
				{
					
				}

				GTree( string model_dir, _int tree_no )
				{
					ifstream fin;
					fin.open( model_dir + "/" + to_string( tree_no ) + ".tree" );

					_int num_nodes;
					fin>>num_nodes;

					for(_int i=0; i<num_nodes; i++)
					{
						GNode* node = new GNode;
						fin>>(*node);
						nodes.push_back(node);
					}
					
					fin.close();
				}

				~GTree()
				{
					for(_int i=0; i<nodes.size(); i++)
						delete nodes[i];
				}

				_float get_ram()
				{
					_float ram = 0;
					for(_int i=0; i<nodes.size(); i++)
						ram += nodes[i]->get_ram();

					return ram;
				}

				void write( string model_dir, _int tree_no )
				{
					ofstream fout;
			        fout.open( model_dir + "/" + to_string( tree_no ) + ".tree" );
					fout<<nodes.size()<<endl;

					for(_int i=0; i<nodes.size(); i++)
					{
						GNode* node = nodes[i];
						fout<<(*node);
					}

					fout.close();
				}
			};

			typedef GTree<Node> Tree;

			Tree* train_tree(SMatF* trn_ft_mat, SMatF* trn_lbl_mat, Param& param, _int tree_no);
			void train_trees( SMatF* trn_X_Xf, SMatF* trn_X_Y, Param& param, string model_dir, _float& train_time );


			SMatF* predict_tree(SMatF* tst_ft_mat, Tree* tree, Param& param);
			SMatF* predict_trees( SMatF* tst_X_Xf, Param& param, string model_dir, _float& prediction_time, _float& model_size );

			_bool optimize_ndcg( SMatF* X_Y, VecI& pos_or_neg );
			void calc_leaf_prob( Node* node, SMatF* X_Y, Param& param );
			_bool optimize_log_loss( SMatF* Xf_X, VecI& y, VecF& C, VecIF& sparse_w, Param& param );
			void setup_thread_locals( _int num_X, _int num_Xf, _int num_Y );
			pairII get_pos_neg_count(VecI& pos_or_neg);
			void test_svm( VecI& X, SMatF* X_Xf, VecIF& w, VecF& values );



// fastXML.cpp
			LOGLVL loglvl = LOGLVL::PROGRESS;  // print progress reports
mutex mtx;	// used to synchronize aggregation of individual tree scores during prediction time
_bool USE_IDCG = true; // if true, optimizes for nDCG; otherwise, optimized for DCG

/* reusable general data containers */ 
thread_local mt19937 reng; // random number generator used during training 
thread_local VecF discounts;
thread_local VecF csum_discounts;
thread_local VecF dense_w;
thread_local VecI countmap;

void setup_thread_locals( _int num_X, _int num_Xf, _int num_Y )
{
    discounts.resize( num_Y );
    csum_discounts.resize( num_Y+1 );
    
    csum_discounts[0] = 1.0;
    _float sumd = 0;
    for( _int i=0; i<num_Y; i++ )
    {
        discounts[i] = 1.0/log2((_float)(i+2));
        sumd += discounts[i];
        
        if(USE_IDCG)
            csum_discounts[i+1] = sumd;
		else
			csum_discounts[i+1] = 1.0;
    }
    dense_w.resize( num_Xf );
    for( _int i=0; i<num_Xf; i++ )
        dense_w[i] = 0;

	countmap.resize( max( num_Xf, num_Y ), 0 );
}

pairII get_pos_neg_count(VecI& pos_or_neg)
{
	pairII counts = make_pair(0,0);
	for(_int i=0; i<pos_or_neg.size(); i++)
	{
		if(pos_or_neg[i]==+1)
			counts.first++;
		else
			counts.second++;
	}
	return counts;
}

typedef signed char schar;
_bool optimize_log_loss( SMatF* Xf_X, VecI& y, VecF& C, VecIF& sparse_w, Param& param )
{
    _int num_X = Xf_X->nr;
    _int num_Xf = Xf_X->nc;
    _int* size = Xf_X->size;
    pairIF** data = Xf_X->data;
    
	_double eps = 0.01;    
	_int l = num_X;
	_int w_size = num_Xf;
	_int newton_iter=0, iter=0;
	_int max_newton_iter = 10;
	_int max_iter = 10;
	_int max_num_linesearch = 20;
	_int active_size;
	_int QP_active_size;

	_double nu = 1e-12;
	_double inner_eps = 1;
	_double sigma = 0.01;
	_double w_norm, w_norm_new;
	_double z, G, H;
	_double Gnorm1_init;
	_double Gmax_old = INF;
	_double Gmax_new, Gnorm1_new;
	_double QP_Gmax_old = INF;
	_double QP_Gmax_new, QP_Gnorm1_new;
	_double delta, negsum_xTd, cond;

    VecD w( num_Xf, 0 );
    VecI index( num_Xf, 0 );
    VecD Hdiag( num_Xf, 0 );
    VecD Grad( num_Xf, 0 );
    VecD wpd( num_Xf, 0 );
    VecD xjneg_sum( num_Xf, 0 );
    VecD xTd( num_X, 0 );
    VecD exp_wTx( num_X, 0 );
    VecD exp_wTx_new( num_X, 0 );
    VecD tau( num_X, 0 );
    VecD D( num_X, 0 );
	
	w_norm = 0;
	for( _int i=0; i<w_size; i++ )
	{
		index[i] = i;

        for( _int j=0; j<size[i]; j++ )
		{
			_int inst = data[i][j].first;
			_float val = data[i][j].second;

			if(y[inst] == -1)
				xjneg_sum[i] += C[inst]*val;
		}
	}

	for( _int i=0; i<l; i++ )
	{
		exp_wTx[i] = exp(exp_wTx[i]);
		_double tau_tmp = 1/(1+exp_wTx[i]);
		tau[i] = C[i]*tau_tmp;
		D[i] = C[i]*exp_wTx[i]*SQ(tau_tmp);
	}

	while(newton_iter < max_newton_iter)
	{
		Gmax_new = 0;
		Gnorm1_new = 0;
		active_size = w_size;

		for(_int s=0; s<active_size; s++)
		{
			_int i = index[s];
			Hdiag[i] = nu;

			_double tmp = 0;
		
			for( _int j=0; j<size[i]; j++ )
			{
				_int inst = data[i][j].first;
				_float val = data[i][j].second;
				Hdiag[i] += SQ(val)*D[inst];
				tmp += val*tau[inst];
			}

			Grad[i] = -tmp + xjneg_sum[i];

			_double Gp = Grad[i]+1;
			_double Gn = Grad[i]-1;
			_double violation = 0;

			if(w[i] == 0)
			{
				if(Gp < 0)
					violation = -Gp;
				else if(Gn > 0)
					violation = Gn;
				//outer-level shrinking
				else if(Gp>Gmax_old/l && Gn<-Gmax_old/l)
				{
					active_size--;
					swap(index[s], index[active_size]);
					s--;
					continue;
				}
			}
			else if(w[i] > 0)
				violation = fabs(Gp);
			else
				violation = fabs(Gn);

			Gmax_new = max(Gmax_new, violation);
			Gnorm1_new += violation;
		}

		if(newton_iter == 0)
			Gnorm1_init = Gnorm1_new;

		if(Gnorm1_new <= eps*Gnorm1_init)
			break;

		iter = 0;
		QP_Gmax_old = INF;
		QP_active_size = active_size;

		for(_int i=0; i<l; i++)
			xTd[i] = 0;

		// optimize QP over wpd
		while(iter < max_iter)
		{
			QP_Gmax_new = 0;
			QP_Gnorm1_new = 0;

			for(_int i=0; i<QP_active_size; i++)
			{
				_llint r = reng();
				_int j = i+r%(QP_active_size-i);
				swap(index[j], index[i]);
			}

			for(_int s=0; s<QP_active_size; s++)
			{
				_int i = index[s];
				H = Hdiag[i];

				G = Grad[i] + (wpd[i]-w[i])*nu;
				for( _int j=0; j<size[i]; j++ )
				{
					_int inst = data[i][j].first;
					_float val = data[i][j].second;
					G += val*D[inst]*xTd[inst];
				}

				_double Gp = G+1;
				_double Gn = G-1;
				_double violation = 0;
				if(wpd[i] == 0)
				{
					if(Gp < 0)
						violation = -Gp;
					else if(Gn > 0)
						violation = Gn;
					//inner-level shrinking
					else if(Gp>QP_Gmax_old/l && Gn<-QP_Gmax_old/l)
					{
						QP_active_size--;
						swap(index[s], index[QP_active_size]);
						s--;
						continue;
					}
				}
				else if(wpd[i] > 0)
					violation = fabs(Gp);
				else
					violation = fabs(Gn);

				QP_Gmax_new = max(QP_Gmax_new, violation);
				QP_Gnorm1_new += violation;

				// obtain solution of one-variable problem
				if(Gp < H*wpd[i])
					z = -Gp/H;
				else if(Gn > H*wpd[i])
					z = -Gn/H;
				else
					z = -wpd[i];

				if(fabs(z) < 1.0e-12)
					continue;
				z = min(max(z,-10.0),10.0);

				wpd[i] += z;

				for( _int j=0; j<size[i]; j++ )
				{
					_int inst = data[i][j].first;
					_float val = data[i][j].second;
					xTd[inst] += val*z;
				}
			}

			iter++;

			if(QP_Gnorm1_new <= inner_eps*Gnorm1_init)
			{
				//inner stopping
				if(QP_active_size == active_size)
					break;
				//active set reactivation
				else
				{
					QP_active_size = active_size;
					QP_Gmax_old = INF;
					continue;
				}
			}

			QP_Gmax_old = QP_Gmax_new;
		}

		delta = 0;
		w_norm_new = 0;
		for(_int i=0; i<w_size; i++)
		{
			delta += Grad[i]*(wpd[i]-w[i]);
			if(wpd[i] != 0)
				w_norm_new += fabs(wpd[i]);
		}
		delta += (w_norm_new-w_norm);

		negsum_xTd = 0;
		for(_int i=0; i<l; i++)
		{
			if(y[i] == -1)
				negsum_xTd += C[i]*xTd[i];
		}

		_int num_linesearch;
		for(num_linesearch=0; num_linesearch < max_num_linesearch; num_linesearch++)
		{
			_double cond = w_norm_new - w_norm + negsum_xTd - sigma*delta;

			for(_int i=0; i<l; i++)
			{
				_double exp_xTd = exp(xTd[i]);
				exp_wTx_new[i] = exp_wTx[i]*exp_xTd;
				cond += C[i]*log((1+exp_wTx_new[i])/(exp_xTd+exp_wTx_new[i]));
			}

			if(cond <= 0)
			{
				w_norm = w_norm_new;
				for(_int i=0; i<w_size; i++)
					w[i] = wpd[i];

				for(_int i=0; i<l; i++)
				{
					exp_wTx[i] = exp_wTx_new[i];
					_double tau_tmp = 1/(1+exp_wTx[i]);
					tau[i] = C[i]*tau_tmp;
					D[i] = C[i]*exp_wTx[i]*SQ(tau_tmp);
				}
				break;
			}
			else
			{
				w_norm_new = 0;
				for(_int i=0; i<w_size; i++)
				{
					wpd[i] = (w[i]+wpd[i])*0.5;

					if(wpd[i] != 0)
						w_norm_new += fabs(wpd[i]);
				}
				delta *= 0.5;
				negsum_xTd *= 0.5;
				for(_int i=0; i<l; i++)
					xTd[i] *= 0.5;
			}
		}

		// Recompute some info due to too many line search steps
		if(num_linesearch >= max_num_linesearch)
		{
			for(_int i=0; i<l; i++)
				exp_wTx[i] = 0;

			for(_int i=0; i<w_size; i++)
			{
				if(w[i]==0) continue;

				for( _int j=0; j<size[i]; j++ )
				{
					_int inst = data[i][j].first;
					_float val = data[i][j].second;
					exp_wTx[inst] += w[i]*val;
				}
			}

			for(_int i=0; i<l; i++)
				exp_wTx[i] = exp(exp_wTx[i]);
		}

		if(iter == 1)
			inner_eps *= 0.25;

		newton_iter++;
		Gmax_old = Gmax_new;
	}

	_float th = 1e-16;
	for( _int i=0; i<w_size; i++ )
	{
		if(fabs(w[i])>th)
			sparse_w.push_back(make_pair(i,w[i]));
        else
            w[i] = 0;
	}
    
    VecF prods( l, 0 );
    for( _int i=0; i<w_size; i++ )
	{
        for( _int j=0; j<size[i]; j++ )
		{
			_int inst = data[i][j].first;
			_float val = data[i][j].second;
            prods[inst] += w[i]*val;
		}
	}
    
    for( _int i=0; i<l; i++ )
        y[i] = prods[i]>=0 ? +1 : -1;
    
	pairII num_pos_neg = get_pos_neg_count(y);

	if(num_pos_neg.first==0 || num_pos_neg.second==0)
	{
		sparse_w.clear();
		return false;
	}

	return true;
}

void calc_leaf_prob( Node* node, SMatF* X_Y, Param& param )
{
    _int lbl_per_leaf = param.lbl_per_leaf;
    _int num_X = X_Y->nc;
    _int num_Y = X_Y->nr;
    _int* size = X_Y->size;
    pairIF** data = X_Y->data;
    
	VecIF& leaf_dist = node->leaf_dist;
	leaf_dist.resize( num_Y );
	for( _int i=0; i<num_Y; i++ )
		leaf_dist[i] = make_pair( i, 0 );

	for( _int i=0; i<num_X; i++ )
	{
		for( _int j=0; j<size[i]; j++ )
        {
			_int lbl = data[i][j].first;
			_float val = data[i][j].second;
			leaf_dist[lbl].second += val;
		}
	}	

	for( _int i=0; i<num_Y; i++ )
		leaf_dist[i].second /= num_X;

	sort( leaf_dist.begin(), leaf_dist.end(), comp_pair_by_second_desc<_int,_float> );
	if( leaf_dist.size()>lbl_per_leaf )
		leaf_dist.resize( lbl_per_leaf );
	sort( leaf_dist.begin(), leaf_dist.end(), comp_pair_by_first<_int,_float> );
}

_bool optimize_ndcg( SMatF* X_Y, VecI& pos_or_neg )
{
    _int num_X = X_Y->nc;
    _int num_Y = X_Y->nr;
    _int* size = X_Y->size;
    pairIF** data = X_Y->data;
    
    _float eps = 1e-6;
    
    VecF idcgs( num_X );
    for( _int i=0; i<num_X; i++ )
        idcgs[i] = 1.0/csum_discounts[ size[i] ];
    
    VecIF pos_sum( num_Y );
    VecIF neg_sum( num_Y );
    VecF diff_vec( num_Y );

    _float ndcg = -2;
    _float new_ndcg = -1;
    
	while(true)
	{
		for(_int i=0; i<num_Y; i++ )
		{
			pos_sum[i] = make_pair(i,0);
			neg_sum[i] = make_pair(i,0);
			diff_vec[i] = 0;
		}

		for( _int i=0; i<num_X; i++ )
		{
			for( _int j=0; j<size[i]; j++ )
			{
				_int lbl = data[i][j].first;
				_float val = data[i][j].second * idcgs[i];

				if(pos_or_neg[i]==+1)
					pos_sum[lbl].second += val;
				else
					neg_sum[lbl].second += val;
			}
		}

		new_ndcg = 0;
		for(_int s=-1; s<=1; s+=2)
		{
			VecIF& sum = s==-1 ? neg_sum : pos_sum;
			sort(sum.begin(), sum.begin()+num_Y, comp_pair_by_second_desc<_int,_float>);

			for(_int i=0; i<num_Y; i++)
			{
				_int lbl = sum[i].first;
				_float val = sum[i].second;
				diff_vec[lbl] += s*discounts[i];
				new_ndcg += discounts[i]*val;
			}
		}
		new_ndcg /= num_X;

		for( _int i=0; i<num_X; i++ )
		{
			_float gain_diff = 0;
            for( _int j=0; j<size[i]; j++ )
			{
                _int lbl = data[i][j].first;
				_float val = data[i][j].second * idcgs[i];
				gain_diff += val*diff_vec[lbl];
			}

			if(gain_diff>0)
				pos_or_neg[i] = +1;
			else if(gain_diff<0)
				pos_or_neg[i] = -1;
		}
	
		if(new_ndcg-ndcg<eps)
			break;
		else
			ndcg = new_ndcg;

	}

	pairII num_pos_neg = get_pos_neg_count(pos_or_neg);
	if(num_pos_neg.first==0 || num_pos_neg.second==0)
		return false;
	return true;
}

void shrink_data_matrices( SMatF* trn_X_Xf, SMatF* trn_X_Y, VecI& n_X, SMatF*& n_trn_Xf_X, SMatF*& n_trn_X_Y, VecI& n_Xf, VecI& n_Y )
{
    trn_X_Xf->shrink_mat( n_X, n_trn_Xf_X, n_Xf, countmap, true ); // countmap is a thread_local variable
    trn_X_Y->shrink_mat( n_X, n_trn_X_Y, n_Y, countmap, false );
}

_bool split_node( Node* node, SMatF* Xf_X, SMatF* X_Y, VecI& pos_or_neg, Param& param )
{
    _int num_X = Xf_X->nr;
	pos_or_neg.resize( num_X );
 
	for( _int i=0; i<num_X; i++ )
	{
		_llint r = reng();

		if(r%2)
			pos_or_neg[i] = 1;
		else
			pos_or_neg[i] = -1;
	}

	// one run of ndcg optimization
	bool success;

	success = optimize_ndcg( X_Y, pos_or_neg );
	if(!success)
		return false;

	VecF C( num_X );
	pairII num_pos_neg = get_pos_neg_count( pos_or_neg );
	_float frac_pos = (_float)num_pos_neg.first/(num_pos_neg.first+num_pos_neg.second);
	_float frac_neg = (_float)num_pos_neg.second/(num_pos_neg.first+num_pos_neg.second);
	_double Cp = param.log_loss_coeff/frac_pos;
	_double Cn = param.log_loss_coeff/frac_neg;  // unequal Cp,Cn improves the balancing in some data sets
	
	for( _int i=0; i<num_X; i++ )
		C[i] = pos_or_neg[i]==+1 ? Cp : Cn;

	// one run of log-loss optimization
	success = optimize_log_loss( Xf_X, pos_or_neg, C, node->w, param );
	if(!success)
		return false;

	return true;
}

void postprocess_node( Node* node, SMatF* trn_X_Xf, SMatF* trn_X_Y, VecI& n_X, VecI& n_Xf, VecI& n_Y )
{
    if( node->is_leaf )
        reindex_VecIF( node->leaf_dist, n_Y );
    else
        reindex_VecIF( node->w, n_Xf );
}

Tree* train_tree( SMatF* trn_X_Xf, SMatF* trn_X_Y, Param& param, _int tree_no )
{
	reng.seed(tree_no);

	_int num_X = trn_X_Xf->nc;
	_int num_Xf = trn_X_Xf->nr;
	_int num_Y = trn_X_Y->nr;

	Tree* tree = new Tree;
	vector<Node*>& nodes = tree->nodes;

	VecI X;
	for(_int i=0; i<num_X; i++)
		X.push_back(i);
	Node* root = new Node( X, 0, param.max_leaf );
	nodes.push_back(root);

	VecI pos_or_neg;

	for(_int i=0; i<nodes.size(); i++)
	{
		if(loglvl == LOGLVL::PROGRESS)
		{
			if(i%1000==0)
				cout<<"\tnode "<<i<<endl;
		}		

		Node* node = nodes[i];
		VecI& n_X = node->X;	
		SMatF* n_trn_Xf_X;
		SMatF* n_trn_X_Y;
		VecI n_Xf;
        VecI n_Y;

		shrink_data_matrices( trn_X_Xf, trn_X_Y, n_X, n_trn_Xf_X, n_trn_X_Y, n_Xf, n_Y );

		if(node->is_leaf)
		{
			calc_leaf_prob( node, n_trn_X_Y, param );
        }
		else
		{
			VecI pos_or_neg;
			bool success = split_node( node, n_trn_Xf_X, n_trn_X_Y, pos_or_neg, param );

			if(success)
			{
				VecI pos_X, neg_X;
				for(_int j=0; j<n_X.size(); j++)
				{
					_int inst = n_X[j];
					if( pos_or_neg[j]==+1 )
						pos_X.push_back(inst);
					else
						neg_X.push_back(inst);
				}
	
				Node* pos_node = new Node( pos_X, node->depth+1, param.max_leaf );
				nodes.push_back(pos_node);
				node->pos_child = nodes.size()-1;

				Node* neg_node = new Node( neg_X, node->depth+1, param.max_leaf );
				nodes.push_back(neg_node);
				node->neg_child = nodes.size()-1;
			}
			else
			{
				node->is_leaf = true;
				i--;
			}
		}
       
        postprocess_node( node, trn_X_Xf, trn_X_Y, n_X, n_Xf, n_Y );

		delete n_trn_Xf_X;
		delete n_trn_X_Y;
	}
	tree->num_Xf = num_Xf;
	tree->num_Y = num_Y;

	return tree;
}

void train_trees_thread( SMatF* trn_X_Xf, SMatF* trn_X_Y, Param param, _int s, _int t, string model_dir, _float* train_time )
{
	Timer timer;
	timer.tic();
    _int num_X = trn_X_Xf->nc;
    _int num_Xf = trn_X_Xf->nr;
    _int num_Y = trn_X_Y->nr;
    setup_thread_locals( num_X, num_Xf, num_Y );
    {
		lock_guard<mutex> lock(mtx);
		*train_time += timer.toc();
    }
    
	for(_int i=s; i<s+t; i++)
	{
		timer.tic();
		cout<<"tree "<<i<<" training started"<<endl;

		Tree* tree = train_tree( trn_X_Xf, trn_X_Y, param, i );
		{
			lock_guard<mutex> lock(mtx);
			*train_time += timer.toc();
		}

		tree->write( model_dir, i );

		timer.tic();
		delete tree;

		cout<<"tree "<<i<<" training completed"<<endl;
		
		{
			lock_guard<mutex> lock(mtx);
			*train_time += timer.toc();
		}
	}
}

void train_trees( SMatF* trn_X_Xf, SMatF* trn_X_Y, Param& param, string model_dir, _float& train_time )
{
	_float* t_time = new _float;
	*t_time = 0;
	Timer timer;
	
	timer.tic();
	trn_X_Xf->append_bias_feat( param.bias );

	_int tree_per_thread = (_int)ceil((_float)param.num_tree/param.num_thread);
	vector<thread> threads;
	_int s = param.start_tree;
	for( _int i=0; i<param.num_thread; i++ )
	{
		if( s < param.start_tree+param.num_tree )
		{
			_int t = min( tree_per_thread, param.start_tree+param.num_tree-s );
			threads.push_back( thread( train_trees_thread, trn_X_Xf, trn_X_Y, param, s, t, model_dir, ref(t_time) ));
			s += t;
		}
	}
	*t_time += timer.toc();	

	for(_int i=0; i<threads.size(); i++)
		threads[i].join();

	train_time = *t_time;
	delete t_time;
}

void test_svm( VecI& X, SMatF* X_Xf, VecIF& w, VecF& values )
{
	values.resize( X.size() );
    _int num_Xf = X_Xf->nr;

	for(_int i=0; i<w.size(); i++)
		dense_w[w[i].first] = w[i].second;

	_int* siz = X_Xf->size;
	pairIF** data = X_Xf->data;

	for(_int i=0; i<X.size(); i++)
	{
		_int inst = X[i];
		_float prod = 0;

		for(_int j=0; j<siz[inst]; j++)
		{
			_int ft = data[inst][j].first;
			_float val = data[inst][j].second;
			prod += val*dense_w[ft];
		}

		values[i] = prod;
	}	

	for(_int i=0; i<w.size(); i++)
		dense_w[w[i].first] = 0;
}

SMatF* predict_tree( SMatF* tst_X_Xf, Tree* tree, Param& param )
{
	_int num_X = tst_X_Xf->nc;
	_int num_Xf = param.num_Xf;
	_int num_Y = param.num_Y;

	vector<Node*>& nodes = tree->nodes;
	Node* node = nodes[0];
	node->X.clear();

	for(_int i=0; i<num_X; i++)
		node->X.push_back(i);

	SMatF* tst_score_mat = new SMatF(num_Y,num_X);
	VecI pos_or_neg( num_X );
	VecF values( num_X );

	for(_int i=0; i<nodes.size(); i++)
	{
		if(loglvl == LOGLVL::PROGRESS)
		{
			if(i%1000==0)
				cout<<"\tnode "<<i<<endl;
		}		

		Node* node = nodes[i];
	
		if(!node->is_leaf)
		{
			VecI& X = node->X;
			test_svm(X, tst_X_Xf, node->w, values);
			for( _int j=0; j<X.size(); j++ )
				pos_or_neg[j] = values[j]>=0 ? +1 : -1;
			Node* pos_node = nodes[node->pos_child];
			pos_node->X.clear();
			Node* neg_node = nodes[node->neg_child];
			neg_node->X.clear();

			for(_int j=0; j<X.size(); j++)
			{
				if(pos_or_neg[j]==+1)
					pos_node->X.push_back(X[j]);
				else
					neg_node->X.push_back(X[j]);
			}
		}
		else
		{
			VecI& X = node->X;
			VecIF& leaf_dist = node->leaf_dist;
			_int* size = tst_score_mat->size;
			pairIF** data = tst_score_mat->data;

			for(_int j=0; j<X.size(); j++)
			{
				_int inst = X[j];
				size[inst] = leaf_dist.size();
				data[inst] = new pairIF[leaf_dist.size()];

				for(_int k=0; k<leaf_dist.size(); k++)
					data[inst][k] = leaf_dist[k];
			}
		}
	}

	return tst_score_mat;
}

void predict_trees_thread( SMatF* tst_X_Xf, SMatF* score_mat, Param param, _int s, _int t, string model_dir, _float* prediction_time, _float* model_size )
{
    Timer timer;
    
    timer.tic();
    _int num_Xf = tst_X_Xf->nr;
    dense_w.resize( num_Xf );
    for( _int i=0; i<num_Xf; i++ )
        dense_w[i] = 0;
	{
		lock_guard<mutex> lock(mtx);
		*prediction_time += timer.toc();
	}

	for(_int i=s; i<s+t; i++)
	{
		if(loglvl == LOGLVL::PROGRESS)
			cout<<"tree "<<i<<" prediction started"<<endl;

		Tree* tree = new Tree( model_dir, i );
        timer.tic();
		SMatF* tree_score_mat = predict_tree( tst_X_Xf, tree, param );

		{
			lock_guard<mutex> lock(mtx);
			score_mat->add(tree_score_mat);
            *model_size += tree->get_ram();
		}

		delete tree;
		delete tree_score_mat;

		if(loglvl == LOGLVL::PROGRESS)
			cout<<"tree "<<i<<" prediction completed"<<endl;
        {
			lock_guard<mutex> lock(mtx);
			*prediction_time += timer.toc();
		}
	}
}

SMatF* predict_trees( SMatF* tst_X_Xf, Param& param, string model_dir, _float& prediction_time, _float& model_size )
{
    _float* p_time = new _float;
	*p_time = 0;

	_float* m_size = new _float;
	*m_size = 0;

	Timer timer;

	timer.tic();
    tst_X_Xf->append_bias_feat( param.bias );
    
    _int num_X = tst_X_Xf->nc;
    SMatF* score_mat = new SMatF( param.num_Y, num_X );

	_int tree_per_thread = (_int)ceil((_float)param.num_tree/param.num_thread);
	vector<thread> threads;

	_int s = param.start_tree;
	for(_int i=0; i<param.num_thread; i++)
	{
		if(s < param.start_tree+param.num_tree)
		{
			_int t = min(tree_per_thread, param.start_tree+param.num_tree-s);
            threads.push_back( thread( predict_trees_thread, tst_X_Xf, ref(score_mat), param, s, t, model_dir, ref( p_time ), ref( m_size ) ));
			s += t;
		}
	}
    *p_time += timer.toc();
	
	for(_int i=0; i<threads.size(); i++)
		threads[i].join();

    timer.tic();

	for(_int i=0; i<score_mat->nc; i++)
		for(_int j=0; j<score_mat->size[i]; j++)
			score_mat->data[i][j].second /= param.num_tree;
    
    model_size = *m_size;
	delete m_size;
    
    *p_time += timer.toc();
	prediction_time = *p_time;
	delete p_time;

	return score_mat;
}



void help()
{
	cerr<<"Sample Usage :"<<endl;
	cerr<<"./fastXML_predict [feature file name] [score file name] [model folder name] T 1 -s 0 -t 50 -q 1"<<endl<<endl;

	cerr<<"-T Number of threads to use. default=[value saved in trained model]"<<endl;
	cerr<<"-s Starting tree index. default=[value saved in trained model]"<<endl;
	cerr<<"-t Number of trees to be grown. default=[value saved in trained model]"<<endl;
	cerr<<"-q quiet option (0/1). default=[value saved in trained model]"<<endl;

	cerr<<"feature and score files are in sparse matrix format"<<endl;
	exit(1);
}

Param parse_param(int argc, char* argv[], string model_folder)
{
	Param param(model_folder+"/param");

	string opt;
	string sval;
	_float val;

	for(_int i=0; i<argc; i+=2)
	{
		opt = string(argv[i]);
		sval = string(argv[i+1]);
		val = stof(sval);

		if(opt=="-T")
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
	check_valid_filename(ft_file,true);
	SMatF* tst_X_Xf = new SMatF(ft_file);

	string score_file = string(argv[2]);
	check_valid_filename(score_file,false);

	string model_folder = string(argv[3]);
	check_valid_foldername(model_folder);

	Param param = parse_param(argc-4, argv+4, model_folder);

	if( param.quiet )
		loglvl = LOGLVL::QUIET;

    _float prediction_time;
    _float model_size;
	SMatF* score_mat = predict_trees( tst_X_Xf, param, model_folder, prediction_time, model_size );

	cout << "prediction time: " << ((prediction_time/tst_X_Xf->nc)*1000.0) << " ms/point" << endl;
	cout << "model size: " << model_size/1e+9 << " GB" << endl;

	score_mat->write(score_file);

	delete tst_X_Xf;
	delete score_mat;
}
