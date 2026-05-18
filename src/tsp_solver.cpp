/*
 Petar 'PetarV' Velickovic
 Linear Programming TSP Solver (Dantzig-Fulkerson-Johnson)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <unistd.h>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <complex>

#include <mst.h>
#include <simplex.h>

#define MAX_N 5001
#define EPS 1e-5
#define PATH_BUFFER 512
#define COMMAND_BUFFER 1024

typedef long long lld;
typedef unsigned long long llu;
typedef unsigned int uint;
using namespace std;

char path_adj[PATH_BUFFER], path_demo[PATH_BUFFER], file_graph[PATH_BUFFER], file_map[PATH_BUFFER];
char path_graph[PATH_BUFFER], path_map[PATH_BUFFER], path_pdf[PATH_BUFFER];
FILE *f_adj, *f_graph, *f_lp;

char cmd[PATH_BUFFER];
char cmd_map[COMMAND_BUFFER];
char cmd_prev[COMMAND_BUFFER];

int n;
double adj[MAX_N][MAX_N];

int simp_n;
double c[MAX_N * MAX_N];

// Constraints
vector<vector<double> > Ap;
vector<double> bp;

int curr_pos = 0;
stack<int> positions;

inline void print_delimiter()
{
    printf("---------------------------------------------------------------\n");
}

inline bool read_token(char *destination, size_t destination_size)
{
    if (destination_size == 0) return false;

    char format[32];
    snprintf(format, sizeof(format), "%%%zus", destination_size - 1);
    return scanf(format, destination) == 1;
}

inline bool read_int(int &value)
{
    return scanf("%d", &value) == 1;
}

inline bool read_two_ints(int &x, int &y)
{
    return scanf("%d%d", &x, &y) == 2;
}

inline bool read_three_ints(int &x, int &y, int &z)
{
    return scanf("%d%d%d", &x, &y, &z) == 3;
}

inline bool join_path(char *destination, size_t destination_size, const char *directory, const char *filename)
{
    const char *separator = "";
    size_t directory_length = strlen(directory);
    if (directory_length > 0 && directory[directory_length - 1] != '/') separator = "/";

    int written = snprintf(destination, destination_size, "%s%s%s", directory, separator, filename);
    return written >= 0 && written < (int)destination_size;
}

inline bool compile_map()
{
    int written = snprintf(cmd_map, sizeof(cmd_map), "cd \"%s\" && pdflatex \"%s\" > /dev/null 2>&1", path_demo, file_map);
    if (written < 0 || written >= (int)sizeof(cmd_map))
    {
        printf("Error: Map compilation command is too long.\n");
        return false;
    }

    int result = system(cmd_map);
    if (result != 0)
    {
        printf("Warning: Map compilation failed. Make sure pdflatex is installed and the demo files are valid.\n");
        return false;
    }

    return true;
}

#ifdef __APPLE__
inline void preview_pdf()
{
    int written = snprintf(cmd_prev, sizeof(cmd_prev), "killall qlmanage > /dev/null 2>&1; qlmanage -p \"%s\" > /dev/null 2>&1 &", path_pdf);
    if (written < 0 || written >= (int)sizeof(cmd_prev))
    {
        printf("Warning: Preview command is too long.\n");
        return;
    }

    printf("%s\n", cmd_prev);

    int result = system(cmd_prev);
    if (result != 0)
    {
        printf("Warning: Preview command failed.\n");
    }
}
#endif

inline int encode_edge(int i, int j)
{
    // Precondition: i > j
    assert(i > j);
    return (((i - 1) * (i - 2)) >> 1) + (j - 1);
}

inline pair<int, int> decode_edge(int x)
{
    int i = 2, j = 1;
    
    while (((i * (i - 1)) >> 1) <= x) i++;
    j = x - (((i - 1) * (i - 2)) >> 1) + 1;
    
    return make_pair(i, j);
}

inline void dump_lp(int n, int m, double **A, double *b, double *c, double v)
{
    fprintf(f_lp, "minimise %.2lf ", v);
    for (int i=0;i<n;i++)
    {
        if (fabs(c[i]) > EPS)
        {
            fprintf(f_lp, "%c %.2lf * x_{%02d, %02d} ", (c[i] >= 0.0 ? '+' : '-'), fabs(c[i]), decode_edge(i).first, decode_edge(i).second);
        }
    }
    fprintf(f_lp, "\nsubject to\n");
    for (int i=0;i<m;i++)
    {
        fprintf(f_lp, "s_{%03d} = %.2lf ", i, b[i]);
        for (int j=0;j<n;j++)
        {
            if (fabs(A[i][j]) > EPS)
            {
                fprintf(f_lp, "%c %.2lf * x_{%02d, %02d} ", (A[i][j] >= 0.0 ? '+' : '-'), fabs(A[i][j]), decode_edge(j).first, decode_edge(j).second);
            }
        }
        fprintf(f_lp, "\n");
    }
}

inline void dump_title(double obj, int n, int m, int iter)
{
    fprintf(f_graph, "\\node[] (title) at (10, 25) {\\Huge Objective value: $%lf$, $%d$ variables, $%d$ constraints, $%d$ iterations};\n", obj, n, m, iter);
}

inline void dump_edge(int x, double val)
{
    fprintf(f_graph, "\\draw[edge,%s] (%d) to node[lab]{%g} (%d);\n", ((fabs(1.0 - val) < EPS) ? "black" : "red"), decode_edge(x).first, val, decode_edge(x).second);
}

inline void write_full_edge(int x, int y)
{
    fprintf(f_graph, "\\draw[edge] (%d) to node[lab]{%d} (%d);\n", x, 1, y);
}

string rep_ext(const char *name, const char *ext)
{
    string ret(name);
    string::size_type ii = ret.find_last_of('.');
    assert(ii != string::npos);
    ret.replace(ii, string::npos, ext);
    
    return ret;
}

int main()
{
    srand(time(NULL));
    
    printf("Linear Programming TSP Solver, implemented by Petar Veličković.\n");
    printf("Thanks to Thomas Sauerwald for the visualisation data!\n");
    print_delimiter();
    
    printf("Enter the path to the file containing the adjacency matrix:\n");
    if (!read_token(path_adj, sizeof(path_adj)))
    {
        printf("Error: Adjacency matrix path could not be read!\n");
        return 1;
    }
    
    if ((f_adj = fopen(path_adj, "r")) == NULL)
    {
        printf("Error: Adjacency matrix file could not be opened!\n");
        return 1;
    }
    
    printf("Processing adjacency matrix...\n");
    if (fscanf(f_adj, "%d", &n) != 1)
    {
        printf("Error: Improper adjacency matrix format!\n");
        return 2;
    }
    
    printf("The graph has %d nodes.\n", n);
    
    for (int i=2;i<=n;i++)
    {
        for (int j=1;j<i;j++)
        {
            if (fscanf(f_adj, "%lf", &adj[i][j]) != 1)
            {
                printf("Error: Improper adjacency matrix format!\n");
                return 3;
            }
        }
    }
    
    fclose(f_adj);
    
    printf("Adjacency matrix successfully read!\n");
    print_delimiter();
    
    printf("Generating basic constraints...\n");
    
    simp_n = (n * (n-1)) >> 1; // the number of edges
    
    // Generating objective function
    for (int i=0;i<simp_n;i++)
    {
        c[i] = -adj[decode_edge(i).first][decode_edge(i).second];
    }
    
    // Generating constraints
    for (int i=1;i<=n;i++)
    {
        vector<double> constr_1(simp_n, 0.0), constr_2(simp_n, 0.0);
        
        for (int j=1;j<=n;j++)
        {
            if (j < i)
            {
                constr_1[encode_edge(i, j)] = -1.0;
                constr_2[encode_edge(i, j)] = 1.0;
            }
            else if (j > i)
            {
                constr_1[encode_edge(j, i)] = -1.0;
                constr_2[encode_edge(j, i)] = 1.0;
            }
        }
        
        Ap.push_back(constr_1);
        bp.push_back(2.0);
        
        Ap.push_back(constr_2);
        bp.push_back(-2.0);
    }
    
    for (int i=0;i<simp_n;i++)
    {
        vector<double> constr(simp_n, 0.0);
        constr[i] = -1.0;
        
        Ap.push_back(constr);
        bp.push_back(1.0);
    }
    
    curr_pos = Ap.size();
    
    printf("Constraints generated!\n");
    print_delimiter();
    
    printf("Enter the path to the folder containing the demo .tex files:\n");
    if (!read_token(path_demo, sizeof(path_demo)))
    {
        printf("Error: Demo folder path could not be read!\n");
        return 4;
    }
    print_delimiter();
    
    printf("Enter the name of the file containing the graph:\n");
    if (!read_token(file_graph, sizeof(file_graph)))
    {
        printf("Error: Graph file name could not be read!\n");
        return 4;
    }
    print_delimiter();
    
    if (!join_path(path_graph, sizeof(path_graph), path_demo, file_graph))
    {
        printf("Error: Graph file path is too long!\n");
        return 4;
    }
    
    printf("Enter the name of the file containing the map:\n");
    if (!read_token(file_map, sizeof(file_map)))
    {
        printf("Error: Map file name could not be read!\n");
        return 4;
    }
    print_delimiter();
    
    if (!join_path(path_map, sizeof(path_map), path_demo, file_map))
    {
        printf("Error: Map file path is too long!\n");
        return 4;
    }
    
    string pdf_file = rep_ext(file_map, ".pdf");
    if (!join_path(path_pdf, sizeof(path_pdf), path_demo, pdf_file.c_str()))
    {
        printf("Error: PDF file path is too long!\n");
        return 4;
    }
    
    
    while (true)
    {
        printf("Enter one of the following:\n");
        printf("- SOLVE : to launch the Simplex Algorithm on the constraints given thus far;\n");
        printf("- REM_LOOP N x_1 x_2 ... x_N : to add a loop-removing constraint for the subset of N nodes x_1, x_2, ..., x_N;\n");
        printf("- REM_LOOP_RNG lo hi : to add a loop-removing constraint for the contiguous interval of nodes [lo, hi];\n");
        printf("- SET x_1 x_2 v : to add constraints that set the edge between x_1 and x_2 to v (where v is either 0 or 1);\n");
        printf("- UNDO N : to remove the N previously generated constraint sets;\n");
        printf("- APPROX_MST : to run the MST-based 2-approximation algorithm;\n");
        printf("- EXIT : to stop the program.\n");
        if (!read_token(cmd, sizeof(cmd)))
        {
            printf("Error: Command could not be read!\n");
            return 5;
        }
        
        if (strcmp(cmd, "SOLVE") == 0)
        {
            int constr_n = curr_pos;
            
            double **A = new double*[constr_n];
            for (int i=0;i<constr_n;i++)
            {
                A[i] = new double[simp_n];
                for (int j=0;j<simp_n;j++)
                {
                    A[i][j] = Ap[i][j];
                }
            }
            
            double *b = new double[constr_n];
            for (int i=0;i<constr_n;i++)
            {
                b[i] = bp[i];
            }
            
            if ((f_lp = fopen("lp.txt", "w")) == NULL)
            {
                printf("Error: LP dump file could not be opened\n");
                return 4;
            }
            
            dump_lp(simp_n, constr_n, A, b, c, 0);
            
            fclose(f_lp);
            
            Simplex *s = new Simplex(simp_n, constr_n, A, b, c, 0);
            auto ret = s -> simplex();
            
            while (std::isnan(get<1>(ret)))
            {
                delete s;
                s = new Simplex(simp_n, constr_n, A, b, c, 0);
                ret = s -> simplex();
            }
            
            delete s;
            for (int i=0;i<constr_n;i++) delete[] A[i];
            delete[] A;
            delete[] b;
            
            printf("Simplex subroutine finished!\n");
            
            vector<double> xs = get<0>(ret);
            double val = get<1>(ret);
            int iters = get<2>(ret);
            
            if (std::isinf(val))
            {
                if (xs[0] == -1) printf("Objective function unbounded!\n");
                else if (xs[0] == -2) printf("Linear program infeasible!\n");
            }
            
            else
            {
                printf("Objective value: %lf\n", val);
                
                //for (int i=0;i<simp_n;i++) if (ret.first[i] != 0) cout << ret.first[i] << endl;
        
                if ((f_graph = fopen(path_graph, "w")) == NULL)
                {
                    printf("Error: Graph file could not be opened!\n");
                    return 4;
                }
                
                dump_title(val, simp_n, constr_n, iters);
                
                for (int i=0;i<simp_n;i++)
                {
                    if (xs[i] > 0) dump_edge(i, xs[i]);
                }
            
                fclose(f_graph);
            
                printf("Results written to the graph file! Compiling the map...\n");
                
                if (compile_map()) printf("Map compiled!\n");
#ifdef __APPLE__
                preview_pdf();
#endif
            }
            
            print_delimiter();
        }
        
        else if (strcmp(cmd, "REM_LOOP") == 0)
        {
            int num;
            if (!read_int(num))
            {
                printf("Error: Loop size could not be read!\n");
                return 5;
            }
            
            vector<int> vals(num);
            for (int i=0;i<num;i++)
            {
                if (!read_int(vals[i]))
                {
                    printf("Error: Loop node could not be read!\n");
                    return 5;
                }
            }
            
            vector<double> constr(simp_n, 0.0);
            
            for (int i=0;i<num;i++)
            {
                for (int j=i+1;j<num;j++)
                {
                    if (vals[i] > vals[j]) constr[encode_edge(vals[i], vals[j])] = -1.0;
                    else if (vals[i] < vals[j]) constr[encode_edge(vals[j], vals[i])] = -1.0;
                }
            }
            
            double val = num - 1.0;
            
            if (curr_pos == (int)Ap.size())
            {
                Ap.push_back(constr);
                bp.push_back(val);
            }
            else
            {
                Ap[curr_pos] = constr;
                bp[curr_pos] = val;
            }
            
            positions.push(curr_pos);
            curr_pos++;
            
            printf("Loop-removing constraint added!\n");
            print_delimiter();
        }
        
        else if (strcmp(cmd, "REM_LOOP_RNG") == 0)
        {
            int lo, hi;
            if (!read_two_ints(lo, hi))
            {
                printf("Error: Loop range could not be read!\n");
                return 5;
            }
            
            assert(hi >= lo);
            
            int num = hi - lo + 1;
            vector<int> vals(num);
            for (int i=0;i<num;i++) vals[i] = lo + i;
            
            vector<double> constr(simp_n, 0.0);
            
            for (int i=0;i<num;i++)
            {
                for (int j=i+1;j<num;j++)
                {
                    if (vals[i] > vals[j]) constr[encode_edge(vals[i], vals[j])] = -1.0;
                    else if (vals[i] < vals[j]) constr[encode_edge(vals[j], vals[i])] = -1.0;
                }
            }
            
            double val = num - 1.0;
            
            if (curr_pos == (int)Ap.size())
            {
                Ap.push_back(constr);
                bp.push_back(val);
            }
            else
            {
                Ap[curr_pos] = constr;
                bp[curr_pos] = val;
            }
            
            positions.push(curr_pos);
            curr_pos++;
            
            printf("Loop-removing constraint added!\n");
            print_delimiter();
        }
        
        else if (strcmp(cmd, "SET") == 0)
        {
            // x(i, j) = v
            int x, y, v;
            if (!read_three_ints(x, y, v))
            {
                printf("Error: Edge setting could not be read!\n");
                return 5;
            }
            
            vector<double> constr1(simp_n, 0.0), constr2(simp_n, 0.0);
            
            if (x > y)
            {
                constr1[encode_edge(x, y)] = -1.0;
                constr2[encode_edge(x, y)] = 1.0;
            }
            else if (x < y)
            {
                constr1[encode_edge(y, x)] = -1.0;
                constr2[encode_edge(y, x)] = 1.0;
            }
            
            positions.push(curr_pos);
            
            if (curr_pos == (int)Ap.size())
            {
                Ap.push_back(constr1);
                bp.push_back(v);
            }
            else
            {
                Ap[curr_pos] = constr1;
                bp[curr_pos] = v;
            }
            
            curr_pos++;
            
            if (v == 1)
            {
                if (curr_pos == (int)Ap.size())
                {
                    Ap.push_back(constr2);
                    bp.push_back(-1.0);
                }
                else
                {
                    Ap[curr_pos] = constr2;
                    bp[curr_pos] = -1.0;
                }
                
                curr_pos++;
            }
            
            printf("Value setting constraints added!\n");
            print_delimiter();
        }
        
        else if (strcmp(cmd, "UNDO") == 0)
        {
            int steps;
            if (!read_int(steps))
            {
                printf("Error: Undo step count could not be read!\n");
                return 5;
            }
            
            bool done = true;
            
            while (steps--)
            {
                if (positions.empty())
                {
                    printf("There's nothing to undo!\n");
                    print_delimiter();
                    done = false;
                    break;
                }
            
                curr_pos = positions.top();
                positions.pop();
            }
            
            if (!done) continue;
            
            printf("Previous generated constraint set(s) removed!\n");
            print_delimiter();
        }
        
        else if (strcmp(cmd, "APPROX_MST") == 0)
        {
            auto mst_sol = tsp_mst(n, adj);
            
            printf("Objective value: %lf\n", mst_sol.first);
            
            //for (int i=0;i<simp_n;i++) if (ret.first[i] != 0) cout << ret.first[i] << endl;
            
            if ((f_graph = fopen(path_graph, "w")) == NULL)
            {
                printf("Error: Graph file could not be opened!\n");
                return 4;
            }
            
            for (uint i=0;i<mst_sol.second.size();i++)
            {
                write_full_edge(mst_sol.second[i].first, mst_sol.second[i].second);
            }
            
            fclose(f_graph);
            
            printf("Results written to the graph file! Compiling the map...\n");
            
            if (compile_map()) printf("Map compiled!\n");
#ifdef __APPLE__
            preview_pdf();
#endif
            print_delimiter();
        }
        
        else if (strcmp(cmd, "EXIT") == 0)
        {
            break;
        }
        
        else
        {
            printf("Incorrectly entered command! Try again.\n");
            print_delimiter();
        }
    }
    
    return 0;
}
