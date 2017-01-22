#include <stdio.h>
#include <iostream>
#include <fstream>
#include <list>
#include <math.h>
#include <cstring>
#include <time.h>

#define INF 100000000    //just infinity
#define _TL 7 // time limit or how many time points track can have
#define _n_c 4 // nubmer of colors for plasmoDB heat_map
#define _d_n 10 // distance between centers of nucleosomes

double Treshold = -1000;
bool TEST_FLAG = false;

using namespace std;

template<class TT> class TMatrix{
      public:
        int m,n;
            TT *p;
            TMatrix(int mm,int nn)
        {
            m=mm;
            n=nn;
            p= new TT[m*n];
        }
        ~TMatrix()
        {
            std::cout<<" Kill TMatrix"<<std::endl;
            delete[] p;
        }
        TT* operator[](int i)
        {
            return  p+n*i;
        }
        void fill(TT x)
        {
            for( int k=0; k < m*n; k++)
                p[ k ] = x;
        }
        void fill(TT x, int k)
        {
            for( int i=0; i < k; i++)
            for( int j = 0; j < k; j++)
                p[ i*n+j ] = x;
        }
        void fillu(TT x, int k)
        {
            for( int i=k; i < m; i++)
            for( int j = k; j < n; j++)
                p[ i*n+j ] = x;
        }
        void Dump(char* str)
        {
            ofstream fout;
            fout.open( str );
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                    fout << p[ n*i + j] << "  ";
                fout << std::endl;
            }
            fout.close();
        }
};

struct Track
{
    public:
        int ids[ _TL ];      // IDs of a assigned nucleosomes
        double pos[ _TL ];  // positon of a track at timepoints
        double prob[ _TL ]; // the score of assignment position-to-track
        double pred;
        double score; //what is the last score on a track
        int cur;  // what is the last timepoint on a track
        int id;   // unique ID of a track
        int time; //time when the track was started
};

inline
void init_labels(TMatrix<double>& cost, double*& lx, double*& ly)
{
    int n = cost.m;
    for (int k=0; k < n; k++)
    {
        lx[k] = 0;
        ly[k] = 0;
    }
    for (int x = 0; x < n; x++)
        for (int y = 0; y < n; y++)
            lx[x] = max(lx[x], cost[x][y]);
}

inline
void add_to_tree(int x, int prevx, TMatrix<double>& cost, double*& lx, double*& ly, bool*& S, int*& prev, double*& slack, int*& slackx)
//x - current vertex,prevx - vertex from X before x in the alternating path,
//so we add edges (prevx, xy[x]), (xy[x], x)
{
    int n = cost.m;
    S[x] = true;                    //add x to S
    prev[x] = prevx;                //we need this when augmenting
    for (int y = 0; y < n; y++)    //update slacks, because we add new vertex to S
        if (lx[x] + ly[y] - cost[x][y] < slack[y])
        {
            slack[y] = lx[x] + ly[y] - cost[x][y];
            slackx[y] = x;
        }
}

inline
void update_labels(int n, double*& slack, double*& lx, double*& ly, bool*& S, bool*& T)
{
    int x, y;
    double  delta = INF;             //init delta as infinity
    for (y = 0; y < n; y++)            //calculate delta using slack
        if (!T[y])
            delta = min(delta, slack[y]);
    for (x = 0; x < n; x++)            //update X labels
        if (S[x]) lx[x] -= delta;
    for (y = 0; y < n; y++)            //update Y labels
        if (T[y]) ly[y] += delta;
    for (y = 0; y < n; y++)            //update slack array
        if (!T[y])
            slack[y] -= delta;
}

void augment(int& max_match, TMatrix<double>& cost, double*& lx, double*& ly, int*& xy, int*& yx)                         //main function of the algorithm
{
    int n = cost.m;
//    std::cout << " Augment " <<100 * max_match / n << "%" <<  std::endl;
    if (max_match == n) return;        //check wether matching is already perfect
    int x, y, root = 0;                    //just counters and root vertex
    int wr = 0, rd = 0;          //q - queue for bfs, wr,rd - write and read
    int* q = new int[n];               //pos in queue
    bool* S = new bool[n];
    bool* T = new bool[n];
    for (int k=0; k < n; k++)
        S[k] = false;       //init set S
    for (int k=0; k < n; k++)//init set T
        T[k] = false;
    int* prev = new int[n];
    for (int k=0; k<n; k++)
        prev[k] = -1;  //init set prev - for the alternating tree
    for (x = 0; x < n; x++)            //finding root of the tree
        if (xy[x] == -1)
        {
            q[wr++] = root = x;
            prev[x] = -2;
	    S[x] = true;
	    break;
        }
    double* slack = new double[n];
    int* slackx = new int[n];
    for (y = 0; y < n; y++)            //initializing slack array
    {
        slack[y] = lx[root] + ly[y] - cost[root][y];
        slackx[y] = root;
    }
    //second part of augment() function
    while (true)                                                        //main cycle
    {
        while (rd < wr)                                                 //building tree with bfs cycle
        {
            x = q[rd++];                                                //current vertex from X part
            for (y = 0; y < n; y++)                                     //iterate through all edges in equality graph
                if (cost[x][y] == lx[x] + ly[y] &&  !T[y])
                {
                    if (yx[y] == -1) break;                             //an exposed vertex in Y found, so
                                                                        //augmenting path exists!
                    T[y] = true;                                        //else just add y to T,
                    q[wr++] = yx[y];                                    //add vertex yx[y], which is matched
                                                                        //with y, to the queue
                    add_to_tree( yx[y], x, cost, lx, ly, S, prev, slack, slackx );//add edges (x,y) and (y,yx[y]) to the tree
                }
            if (y < n) break;                                           //augmenting path found!
        }
        if (y < n) break;                                               //augmenting path found!

        update_labels(n, slack, lx, ly, S, T);                          //augmenting path not found, so improve labeling
        wr = rd = 0;
        for (y = 0; y < n; y++)
        //in this cycle we add edges that were added to the equality graph as a
        //result of improving the labeling, we add edge (slackx[y], y) to the tree if
        //and only if !T[y] &&  slack[y] == 0, also with this edge we add another one
        //(y, yx[y]) or augment the matching, if y was exposed
            if (!T[y] &&  slack[y] == 0)
            {
                if (yx[y] == -1)                                        //exposed vertex in Y found - augmenting path exists!
                {
                    x = slackx[y];
                    break;
                }
                else
                {
                    T[y] = true;                                        //else just add y to T,
                    if (!S[yx[y]])
                    {
                        q[wr++] = yx[y];                                //add vertex yx[y], which is matched with
                                                                        //y, to the queue
                        add_to_tree(yx[y], slackx[y], cost, lx, ly, S, prev, slack, slackx); //and add edges (x,y) and (y,
                                                                        //yx[y]) to the tree
                    }
                }
            }
        if (y < n) break;                                               //augmenting path found!
    }
    delete[] q;
    delete[] T;
    delete[] S;
    delete[] slack;
    delete[] slackx;
    if (y < n)                                                          //we found augmenting path!
    {
        max_match++;                                                    //increment matching
        //in this cycle we inverse edges along augmenting path
        for (int cx = x, cy = y, ty; cx != -2; cx = prev[cx], cy = ty)
        {
            ty = xy[cx];
            yx[cy] = cx;
            xy[cx] = cy;
        }
        augment(max_match, cost, lx, ly, xy, yx);                                                      //recall function, go to step 1 of the algorithm
    }
    delete[] prev;
}


struct Nucleosome
{
    public:
        double pos;
        double size;
        double score;
        double scoreL;
        double scoreR;
        double var;
        int id;
};

inline
bool compare_pos (Nucleosome x, Nucleosome y)
{
    if (x.pos >= y.pos) return false;
    return true;
}

inline
bool compare_pred (Track x, Track y)
{
    if (x.pred >= y.pred) return false;
    return true;
}

inline
bool same_pos (Nucleosome x, Nucleosome y)
{
    if ( fabs(x.pos - y.pos) < _d_n ) return true;
    return false;
}

inline
double Print ( std::list<Nucleosome> L, int Id, double& size, double& var, double& score, double& scoreL, double& scoreR )
{
    std::list<Nucleosome>::iterator it;
    if ( Id >= L.size() )
    {
        return -1;
        size = -1;
        var = -1;
        score = -1;
        scoreL = -1;
        scoreR = -1;
    }
    for( it = L.begin(); it != L.end(); it++)
    {
        if ( it->id == Id ) break;
    }
    size = it->size;
    var = it->var;
    score = it->score;
    scoreL = it->scoreL;
    scoreR = it->scoreR;
    return it->pos;
}

double nuc_pos ( std::list<Nucleosome> L, int Id )
{
    std::list<Nucleosome>::iterator it;
    if ( Id >= L.size() || Id < 0) return -1;
    for( it = L.begin(); it != L.end(); it++)
    {
        if ( it->id == Id ) break;
    }
    return it->pos;
}


double nuc_score ( std::list<Nucleosome> L, int Id )
{
    std::list<Nucleosome>::iterator it;
    if ( Id >= L.size() || Id < 0) return 0;
    for( it = L.begin(); it != L.end(); it++)
    {
        if ( it->id == Id ) break;
    }
    return it->score;
}

void Load_NucleosomeTF ( std::list<Nucleosome>& L, char* str ) //Template Filtering output loading
{
    FILE  *f_in;
    f_in=fopen( str , "r");
    double x, y, z, z1, var, score, scoreL, scoreR, varL, varR;
    int id_count=0;
    Nucleosome T;
    while ( !feof(f_in) )
    {
        fscanf( f_in, "%lf %lf %lf ", &x,  &y, &score);
//            z = ( x + y ) / 2;
//        z1 = fabs( y - x );
        T.pos = x;
        T.var = 0;
        T.size = y;
        T.score = score;
        T.scoreL = 0;
        T.scoreR = 0;
        T.id = id_count;
        if ( y > 40 )
        {
            id_count++;
            L.push_back( T );
//        std::cout<< x << " " << var << " " << y << " " << score << " " << scoreL << " " << scoreR << std::endl;
        }
//        std::cout<< id_count << std::endl;
//      if ( z1 != 146) cout<< z <<" "<<z1<<endl;
    }
//    L.sort(compare_pos);
//    L.unique(same_pos);
    fclose( f_in );
    std::cout << str << " loaded; " << std::endl;
}

void Load_Tracks2 ( std::list<Nucleosome>& L, char* str, int column ) //Template Filtering output loading
{
    FILE  *f_in;
    f_in=fopen( str , "r");
    double x[2], y[2];
    int id_count=0;
    Nucleosome T;
    std::cout << "start scanning...";
    while ( !feof(f_in) )
    {
        fscanf( f_in, "%lf/%lf %lf/%lf", &x[0], &y[0], &x[1], &y[1]);
//            z = ( x + y ) / 2;
//        z1 = fabs( y - x );
        // std::cout << x[0] << x[1] << y[0] << y[1];
        T.pos = x[column];
        T.var = 0;
        T.size = 1;
        T.score = 1;
        T.scoreL = 0;
        T.scoreR = 0;
        T.id = id_count;
        if ( y[column] > 0 )
        {
            id_count++;
            L.push_back( T );
//        std::cout<< x << " " << var << " " << y << " " << score << " " << scoreL << " " << scoreR << std::endl;
        }
//        std::cout<< id_count << std::endl;
//      if ( z1 != 146) cout<< z <<" "<<z1<<endl;
    }
//    L.sort(compare_pos);
//    L.unique(same_pos);
    fclose( f_in );
    std::cout << str << " loaded; " << std::endl;
}

void Load_Tracks3 ( std::list<Nucleosome>& L, char* str, int column ) //Template Filtering output loading
{
    FILE  *f_in;
    f_in=fopen( str , "r");
    double x[3], y[3];
    int id_count=0;
    Nucleosome T;
    std::cout << "start scanning...";
    while ( !feof(f_in) )
    {
        fscanf( f_in, "%lf/%lf %lf/%lf %lf/%lf", &x[0], &y[0], &x[1], &y[1], &x[2], &y[2]);
//            z = ( x + y ) / 2;
//        z1 = fabs( y - x );
        // std::cout << x[0] << x[1] << y[0] << y[1];
        T.pos = x[column];
        T.var = 0;
        T.size = 1;
        T.score = 1;
        T.scoreL = 0;
        T.scoreR = 0;
        T.id = id_count;
        if ( y[column] > 0 )
        {
            id_count++;
            L.push_back( T );
//        std::cout<< x << " " << var << " " << y << " " << score << " " << scoreL << " " << scoreR << std::endl;
        }
//        std::cout<< id_count << std::endl;
//      if ( z1 != 146) cout<< z <<" "<<z1<<endl;
    }
//    L.sort(compare_pos);
//    L.unique(same_pos);
    fclose( f_in );
    std::cout << str << " loaded; " << std::endl;
}

void Load_Tracks4 ( std::list<Nucleosome>& L, char* str, int column ) //Template Filtering output loading
{
    FILE  *f_in;
    f_in=fopen( str , "r");
    double x[4], y[4];
    int id_count=0;
    Nucleosome T;
    std::cout << "start scanning...";
    while ( !feof(f_in) )
    {
        fscanf( f_in, "%lf/%lf %lf/%lf %lf/%lf %lf/%lf", &x[0], &y[0], &x[1], &y[1], &x[2], &y[2], &x[3], &y[3]);
//            z = ( x + y ) / 2;
//        z1 = fabs( y - x );
        // std::cout << x[0] << x[1] << y[0] << y[1];
        T.pos = x[column];
        T.var = 0;
        T.size = 1;
        T.score = 1;
        T.scoreL = 0;
        T.scoreR = 0;
        T.id = id_count;
        if ( y[column] > 0 )
        {
            id_count++;
            L.push_back( T );
//        std::cout<< x << " " << var << " " << y << " " << score << " " << scoreL << " " << scoreR << std::endl;
        }
//        std::cout<< id_count << std::endl;
//      if ( z1 != 146) cout<< z <<" "<<z1<<endl;
    }
//    L.sort(compare_pos);
//    L.unique(same_pos);
    fclose( f_in );
    std::cout << str << " loaded; " << std::endl;
}

void Load_Tracks6 ( std::list<Nucleosome>& L, char* str, int column ) //Template Filtering output loading
{
    FILE  *f_in;
    f_in=fopen( str , "r");
    double x[6], y[6];
    int id_count=0;
    Nucleosome T;
    std::cout << "start scanning...";
    while ( !feof(f_in) )
    {
        fscanf( f_in, "%lf/%lf %lf/%lf %lf/%lf %lf/%lf %lf/%lf %lf/%lf", &x[0], &y[0], &x[1], &y[1], &x[2], &y[2], &x[3], &y[3], &x[4], &y[4],  &x[5], &y[5]);
//            z = ( x + y ) / 2;
//        z1 = fabs( y - x );
        // std::cout << x[0] << x[1] << y[0] << y[1];
        T.pos = x[column];
        T.var = 0;
        T.size = 1;
        T.score = 1;
        T.scoreL = 0;
        T.scoreR = 0;
        T.id = id_count;
        if ( y[column] > 0 )
        {
            id_count++;
            L.push_back( T );
//        std::cout<< x << " " << var << " " << y << " " << score << " " << scoreL << " " << scoreR << std::endl;
        }
//        std::cout<< id_count << std::endl;
//      if ( z1 != 146) cout<< z <<" "<<z1<<endl;
    }
//    L.sort(compare_pos);
//    L.unique(same_pos);
    fclose( f_in );
    std::cout << str << " loaded; " << std::endl;
}




void Load_Nucleosome ( std::list<Nucleosome>& L, char* str )
{
    FILE  *f_in;
    f_in=fopen( str , "r");
    double x, y, z, z1, var, score, scoreL, scoreR, varL, varR;
    int id_count=0;
    Nucleosome T;
    while ( !feof(f_in) )
    {
        fscanf( f_in, "%lf %lf %lf %lf %lf %lf %lf %lf", &x, &var, &varL, &varR, &y, &score, &scoreL, &scoreR);
//            z = ( x + y ) / 2;
//        z1 = fabs( y - x );
        T.pos = x;
        T.var = var;
        T.size = y;
        T.score = score;
        T.scoreL = scoreL;
        T.scoreR = scoreR;
        T.id = id_count;
        if ( y > 30 )
        {
            id_count++;
            L.push_back( T );
//        std::cout<< x << " " << var << " " << y << " " << score << " " << scoreL << " " << scoreR << std::endl;
        }
//        std::cout<< id_count << std::endl;
//      if ( z1 != 146) cout<< z <<" "<<z1<<endl;
    }
//    L.sort(compare_pos);
//    L.unique(same_pos);
    fclose( f_in );
    std::cout << str << " loaded; " << std::endl;
}

inline
void update_pred(std::list<Track>::iterator iter)
{
    double score = 0;
    double res = 0;
    int count = -1;
	for( int i = 0; i <= iter->cur ; i++)
    {
        if ( iter->pos[ i ] > 0)
        {
            res = iter -> pos[ i ];
            count++;
        }
    }
    if (count >= 0)
    {
        iter->pred = res;
    }
    else
    {
        iter->pred = -1;
        iter->score = 0;
    }
}

inline
void Predict( std::list<Track> T)
{
	std::list<Track>::iterator it;
	for( it = T.begin(); it != T.end(); it++)
		update_pred(it);
}

inline
void Fill_Cost_Dummy( TMatrix<double>& Res, std::list<Track> T, std::list<Nucleosome> N, int boundary)
{
//    int tsize = T.size();
//    int nsize = N.size();
//    int size = Res.m;
    std::list<Track>::iterator it;
    for (it = T.begin(); it != T.end(); it++)
        {
            Res[ it->id ][ it->id + boundary ] = (it->score+1) ;
        }

    std::list<Nucleosome>::iterator jt;
    for (jt = N.begin(); jt != N.end(); jt++)
        {
            Res[ jt->id + boundary ][ jt->id ] = (jt->score+1) ;
        }
}

inline
void Fill_Cost( TMatrix<double>& Res, std::list<Track> T, std::list<Nucleosome> N, double tresh)
{
//    int tsize = T.size();
//    int nsize = N.size();
//    int size = Res.m;
    std::list<Track>::iterator it;
    std::list<Nucleosome>::iterator jt;
    for (it = T.begin(); it != T.end(); it++)
        for (jt = N.begin(); jt != N.end(); jt++)
        {
           Res[ it->id ][ jt->id ] = - ( it->pred - jt->pos ) * ( it->pred - jt->pos );   // Distance function
    //        Res[ it->id ][ jt->id ] = 3000/( 1000+(it->pred - jt->pos)*(it->pred - jt->pos)) * (jt->score+1) * (it->score + 1);
        }

}

void Track_to_PlasmoDB(std::list<Track> T, char* str, char* chr)
{
    std::list<Track>::iterator it;
    ofstream fout;
    fout.open( str );
//    fout<< "reference = Pf3D7_" << chr << std::endl;
    for( int j=1; j<= _TL; j++)
    {
        fout<< "[TIME_" << 6*(j-1) << "]" << std::endl;
        fout<< "glyph = heat_map" << std::endl;
        fout<< "start_color = blue" << std::endl;
        fout<< "end_color = red" << std::endl;
        fout<< "min_score = 0" << std::endl;
        fout<< "max_score = " << _n_c - 1 << std::endl;
    }
/*    fout<< "[PRED]" << std::endl;
    fout<< "glyph = heat_map" << std::endl;
    fout<< "start_color = green" << std::endl;
    fout<< "end_color = yellow" << std::endl;
    fout<< "min_score = 0" << std::endl;
    fout<< "max_score = " << _n_c - 1 << std::endl;
*/
    int cc = 0;
    for( it = T.begin(); it != T.end(); it++)
    {
        TEST_FLAG=false;
        if (it->id == 3205)
        {
            TEST_FLAG = true;
        }
        cc++;
        for( int i=0; i <= it->cur; i++)
        {
            if ( it->pos[i] >  0 )
            {
            fout<< "TIME_" << 6*( i + it->time ) << " ";
            fout<< "\"" << it->id << "\" " << " Pf3D7_" << chr << ":" << (int)it->pos[i]-60 << ".." << (int)it->pos[i]+60 << " ";
            fout<< "score=" << cc % _n_c << std::endl;
            }
        }
//        fout<< "PRED ";
//        fout<< "\"" << it->id << "\"" << " Pf3D7_" << chr << ":" << (int)it->pred-73 << ".." << (int)it->pred+73 << " ";
//        fout<< "score=" << it->id % _n_c << std::endl;
    }
    fout.close();
}


void Track_Nucleosome_Print(std::list<Track> T, char* str)
{
    std::list<Track>::iterator it;
    ofstream fout;
    fout.open( str );
    for( it = T.begin(); it != T.end(); it++)
    {
        fout<< it->id << " ";
        for(int i = 0; i < it->time; i++)
            fout << " -1 ";
        for(int i = 0; i <= it->cur; i++)
        {
            fout <<" "<< it->pos[ i ];
        }
        fout << std::endl;
    }
    fout.close();
}

void Track_Print(std::list<Track> T, char* str)
{
    std::list<Track>::iterator it;
    ofstream fout;
    fout.open( str );
    int perfect_hit = 0;
    double perfect_score = 0;
    int miss_hit = 0;
    for( it = T.begin(); it != T.end(); it++)
    {
        // fout<< it->id << " * ";
        for(int i = 0; i < it->time; i++)
            fout << " -1 ";
        for(int i = 0; i <= it->cur; i++)
        {
            if ( it->pos[ i ] > 0 && i > 0)
            {
                perfect_hit++;
                perfect_score += it->prob[ i ];
            }
            else
            {
                miss_hit++;
            }
            fout <<" "<< it->pos[i] <<" " ;

        }
        fout << std::endl;
        miss_hit += _TL - it->cur - 1;
    }
//    fout << "Perfect match: " << perfect_hit << std::endl << "Score: " << ( perfect_score + 0.0 ) / perfect_hit << std::endl;
//    fout << "Missed hit: " << miss_hit;
    fout.close();
}

inline
void Track_Prolong(std::list<Track>& Tracks, std::list<Nucleosome>& N1, std::list<Nucleosome>& N2, int*& ass, TMatrix<double>& Cost)
{
    std::list<Track>::iterator T_it;
    std::list<Nucleosome>::iterator N_it;
    Track Tr;
    int as, x, c, t = 0;
    double p;
    for( T_it = Tracks.begin(); T_it != Tracks.end(); T_it++)
    {
        c = T_it->cur + 1;    // new curent lenght of a track
        if ( c > t)
            t = c;
        as = ass[ T_it->id ]; // id of a new nucleosome in a track
        x = T_it->id;         // id of a track
        p = nuc_pos( N2, as); // position of a new nucleosome in a track
        double score_new = nuc_score( N2, as);
        T_it->cur = c;        // update the lenth of a track
        T_it->prob[ c ] = Cost[ x ][ as ]; // strore score of a new assignment
        T_it->ids[ c ] = as;  // store the id of added nucleosome
        T_it->pos[ c ] = p;   // store position of added nucleosome
        T_it->score = score_new; //store new track score
        if ( p > 0)
            T_it->pred = p;       // NB! the virtual track position

        update_pred(T_it);
    }

    int tr_size = Tracks.size();
    int pr_size = Cost.m;
    int id_count = tr_size;
    for( int i = tr_size; i < pr_size; i++)
    {
        as = ass[ i ];
        p = nuc_pos( N2, as);
        if ( ( as > 0 ) && ( p > 0 ))
        {
            Tr.id = id_count;
            id_count++;
            Tr.cur = c = 0;
            Tr.ids[ c ] = as;
            Tr.pos[ c ] = p;
            Tr.prob[ c ] = Cost[ i ][ as ];
            Tr.time = t;
            Tr.pred = p;      // NB!
            Tracks.push_back( Tr );
        }
    }

}

void Assignment_Problem( std::list<Track>& Tr, std::list<Nucleosome> Nu1, std::list<Nucleosome> Nu2)
{
    int Number_bad = max ( Tr.size(), Nu2.size() ) ;
    int problem_size = Number_bad +  max ( Tr.size(), Nu2.size() );
    TMatrix<double> Cost( problem_size, problem_size);
    Cost.fill( Treshold );
    Cost.fill( 0 , problem_size - Number_bad);
    Cost.fillu( 0, problem_size - Number_bad );
    Fill_Cost( Cost, Tr, Nu2, Treshold);
   // Fill_Cost_Dummy( Cost, Tr, Nu2, problem_size - Number_bad);
    //Cost.Dump("CostDump.txt");
    std::cout << " Cost matrix is ready " << std::endl;
    double* lx = new double[problem_size];
    double* ly = new double[problem_size];
    int max_match = 0;                    //number of vertices in current matching
    int* xy = new int[problem_size];
    int* yx = new int[problem_size];
    for (int k=0 ; k < problem_size; k++)
    {
    xy[k] = -1;
    yx[k] = -1;
    }
    init_labels(Cost, lx, ly);                    //step 0
    augment(max_match, Cost, lx,  ly, xy, yx);   //steps 1-3
    // ofstream f_out;
    // f_out.open( str );
    double ret = 0;                      //weight of the optimal matching
    for (int x = 0; x < problem_size ; x++)       //forming answer there
    {
        ret += Cost[x][xy[x]];
//        f_out << "( " << Print( Nu1, x) << " , " << Print( Nu2, xy[x] ) << " ) " << Cost[ x ][ xy[x] ] << std::endl ;
        double size, var, score, scoreL, scoreR;
        // f_out << "( " << x << " " << Print( Nu1, x, size, var, score, scoreL, scoreR) << " , " << xy[x] << " " << Print( Nu2, xy[x] , size, var, score, scoreL, scoreR) << " ) " << Cost[ x ][ xy[x] ] << std::endl ;
    }
    // f_out << ret << std::endl;
    cout << ret << std::endl;
//    Track_Print( Tr, "track_dump1.txt");
    Track_Prolong( Tr, Nu1, Nu2, xy, Cost);
//    Track_Print( Tr, "track_dump2.txt");
    Tr.sort( compare_pred );
    if ( Tr.size() == problem_size )
        std::cout<< "Limit for dummy nodes exceeded" << std::endl;
    std::cout<< "Track update done" << std::endl;

    delete[] lx;
    delete[] ly;
    delete[] xy;
    delete[] yx;

}

int main(int argc, char *argv[])
{
    time_t before,after ;
    before = time( NULL );
    std::cout<< before<<std::endl;
    list<Nucleosome> N_0, N_1, N_2, N_3, N_4, N_5, N_6, N_7;
    Load_Tracks6 (N_0, argv[1], 0 );
    Load_Tracks6 (N_1, argv[1], 1 );
    Load_Tracks6 (N_2, argv[1], 2 );
    Load_Tracks6 (N_3, argv[1], 3 );
    Load_Tracks6 (N_4, argv[1], 4 );
    Load_Tracks6 (N_5, argv[1], 5 );
    N_0.sort( compare_pos );
    std::cout << "size of list 1: " << N_0.size() << std::endl;
    N_1.sort( compare_pos);
    std::cout << "size of list 2: " << N_1.size() << std::endl;
    N_2.sort( compare_pos);
    std::cout << "size of list 3: " << N_2.size() << std::endl;
    N_3.sort( compare_pos);
    std::cout << "size of list 3: " << N_3.size() << std::endl;
    N_4.sort( compare_pos);
    std::cout << "size of list 3: " << N_4.size() << std::endl;
    std::list<Track> Tracks;
    Track Tr;
    std::list<Nucleosome>::iterator it, jt;
    it = N_0.begin();
    int count = 0;
    while(it != N_0.end())
    {
        Tr.ids[ 0 ] = it->id;
        Tr.cur = 0;
        Tr.pred = it->pos;
        Tr.pos[ 0 ] = it->pos;
        Tr.id = count;
        Tr.prob[0] = 0;
        Tr.time = 0;
        Tracks.push_back(Tr);
        count++;
        it++;
    }
//    Track_Print( Tracks, "track_in.txt");
    Assignment_Problem(Tracks, N_0, N_1 );
    Assignment_Problem(Tracks, N_0, N_2 );
    Assignment_Problem(Tracks, N_0, N_3 );
    Assignment_Problem(Tracks, N_0, N_4 );
    Assignment_Problem(Tracks, N_0, N_5 );
//    Track_to_PlasmoDB(Tracks, "plasmo1.txt", argv[12]);
//    Track_Print( Tracks, "OUTPUT1/track1.txt" );

    Track_Print( Tracks, argv[2]);
//    Track_Nucleosome_Print(Tracks, "track_debug");
    after = time( NULL );
    std::cout<<after<<std::endl;
    std::cout<<" took time : "<<after-before<<std::endl;
    return 0;
}
