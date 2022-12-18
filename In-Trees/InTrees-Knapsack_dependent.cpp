#include <bits/stdc++.h>
#include <ctime>
#include<math.h>
#define ll int64_t
#define MX 1000007
#define MX1 1003
using namespace std;

// The path for the KP instances
char PATH[100]="E:\\PaperCOR\\1.Experim\\0.Inst\\1.C_Uncor\\KP_";//

//coefficient to compute the number of nodes, and the sum
float coeff[5]={0.1,0.2,0.5,0.8};


//----------------------- print the elements of array, vector and vector of vectors ---------
void printVect(vector<int> arr) 
{ 
    cout<<"\n";
	for (int i = 0; i < arr.size(); i++) 
        cout << arr[i] << " "; 
}
void printVectVect(vector<vector<int> > K)
{
	cout<<"\n-----------\n";
		for(int i=0;i<K.size();i++)
		{
			cout<<i<<"--> ";
			for(int j=0;j<K[i].size();j++)
			{
				cout<<K[i][j]<<" ";
			}
			cout<<"\n";
		}
}
//------------------------- affiche the items filled from the file -----------------------
void display_fillin (vector<vector<int> > weights, vector<vector<int> > profits, int m)
{
	cout<<"\n distribution of items weights and profits:";
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<weights[i].size();j++)
			cout<<weights[i][j]<<", "<<profits[i][j]<<"\n";;
		cout<<"\n";	
	}
}
	
//----------------------- print the elements of an array ---------
void printArr(int arr[], int n) 
{ 
    cout<<"\n";
	for (int i = 0; i < n; i++) 
        cout << arr[i] << " "; 
}
void printMat(vector<vector<int> > K,int n, int W)
{
	cout<<"\n";
	for(int i=0;i<=n;i++)
	{
    	for(int j=0;j<=W;j++)
      		cout<<"  "<<K[i][j];
    	cout<<endl;
	}
}

//------------ print prufer squence ------------ 
void print_fruper(int s[], int n) 
{ 
	cout<<"\n Prufer squence: ";
	cout<<"( ";
    for (int i = 0; i < n; i++) 
    {
    	cout << s[i];
    	if(i!=n-1)
			cout <<  ", "; 
	}
	cout<<") \n";
} 
//-------------- display the subtrees of the tree------------------
void display_subtrees(vector<vector<int> > subtrees, int root)
{
	cout<<"\nRoot: "<<root<<"\n";
	cout<<"Nbr_subtrees = "<<subtrees.size()<<"\n";
	for (int i = 0; i < subtrees.size(); i++)
    {
        for (int j = 0; j < subtrees[i].size(); j++)
        {
            cout << subtrees[i][j] << " ";
        }    
        cout << endl;
    }
}
//-------------- display the subtrees of the tree without root------------------
void display_subtrees(vector<vector<int> > subtrees)
{
	cout<<"\nNbr_subtrees = "<<subtrees.size()<<"\n";
	for (int i = 0; i < subtrees.size(); i++)
    {
        for (int j = 0; j < subtrees[i].size(); j++)
        {
            cout << subtrees[i][j] << " ";
        }    
        cout << endl;
    }
}

//-------------------generate a list of number which are coefficient of another number
void list_nbr(int l[],float coeff[],int size,int n)
{
	for(int i=0;i<size;i++)
		l[i]=floor(n*coeff[i]); 
}
//-------------------- generates a random number in [a,b] ---------------
int rand_a_b(int a,int b) 
{
    int r;
    b; //include b
    if(a>=b)
        r = a;
    else
        r = rand()%(b-a) + a;
    return r;
}
// ---------- Compute the max value in a vector v -----------------------------------------
double max_v(vector<double> x) 
{
    double m = x[0];
    for(int i=1;i<x.size();i++)
        if(x[i] > m)
            m = x[i];
    return m;
}
//----------- Compute the average of the values in a vector v----------
double avg_v(vector<double> x) 
{
    double s,m;
    s=0;
    for(int i=0;i<x.size();i++)
        s = s+ x[i];

    m = s/2;//m = s/x.size();
    return m;
}
//---------------------------------- Max of two integer-------------------------------
int max(int x, int y) {
   return (x > y) ? x : y;
}
//---------------------- generate a list of m random non-negative integers whose sum is n 
void randomList(int m, int n,int arr[]) 
{ 
    //initialze the m number to 1 to ensure that the random are all positive (!=0)
	for (int i = 0; i < m; i++)
    	arr[i]=1;
    	
	srand(time(0)); 
       
	// To make the sum of the final list as n 
    for (int i = 0; i < n-m; i++) { 
  
        // Increment any random element from the array by 1
        arr[rand() % m]++; 
    } 
} 

//------------- Get the file name ----------------
void f_name(char file_name[100],int n, int inst)
{
	char ch[20];
		
	// concatenate KP and n
    //strcpy(file_name,PATH);strcat(file_name,"KP_"); 
	strcpy(file_name,PATH);
	sprintf(ch,"%d",n);
	strcat(file_name,ch);
	
	// concatenate KP_n and number of instance 
    strcat(file_name,"_");
	sprintf(ch,"%d",inst);
	strcat(file_name,ch);
	
	// Add extension "txt"
    strcat(file_name,".txt");
}

//-------------------- Load data from file -----------------------------
void load_data_file(char file_name[50],int* t_W,vector<vector<int> > &weights, vector<vector<int> > &profits, int arr[]) 
{
    int i,j;
    FILE* fichier = NULL;
    int x,W;
    
    char ch[20];
	fichier = fopen(file_name, "r");
	
	if(fichier!= NULL)
    {
        fscanf(fichier, "%d", &x); // the number of items
        fscanf(fichier, "%d", &W); // the total weight of all items
        fscanf(fichier, "%d", &x);
    }

	if(fichier!= NULL)
     {
     	printf("\n Open successes %s",file_name);
        i = 0;
        j = 0;
    	do
		{
			if(j < arr[i])
            {
                weights[i].push_back(x); 
				fscanf(fichier, "%d", &x);
				
				profits[i].push_back(x);
				fscanf(fichier, "%d", &x);
            }
            else
            {
                i++;
                j = -1;
            }
            j++;
            
        }while(!feof(fichier));
        fclose(fichier);
    }
    else
    {
        printf("\n Impossible to open the file %s",file_name);
    }
    *t_W=W;
}

//---------------generate a random Prufer sequence of n-2 integers-----------
void fruper(int n, int l, int s[]) //the corresponding tree will have n nodes and l levels
{ 
	int i;
	srand(time(0)); 

    // If the level is not specified 
	if(l ==0)
	{
		// generate a list of n-2 random integer
		for (i = 0; i <n-2; i++)
			s[i] = 1+(rand() % (n-1));
	}
	else
	{
		// ensure distinct values for the first l cases, then generate random values (0...l) for the rest cases
		for (i = 0; i <n-2; i++)
			if(i<l)
				s[i] = i+1;
			else
				s[i] = 1+ (rand() % (l-1));
	}
} 

//--------------------- Generate edges of tree represented by a given Prufer code -----------
void printTreeEdges(int prufer[], int m, int edges[][2]) 
{ 
    int l=0;
	int vertices = m + 2; 
    int vertex_set[vertices]; 
    //int edges[vertices-1][2];
  
    // Initialize the array of vertices 
    for (int i = 0; i < vertices; i++) 
        vertex_set[i] = 0; 
  
    // Number of occurrences of vertex in code 
    for (int i = 0; i < vertices - 2; i++) 
        vertex_set[prufer[i] - 1] += 1; 
  
    int j = 0; 
    for (int i = 0; i < vertices - 2; i++) { 
        for (j = 0; j < vertices; j++) { 
            // If j+1 is not present in prufer set 
            if (vertex_set[j] == 0) { 
                // Remove from Prufer set and print pair. 
                vertex_set[j] = -1; 
                //cout << "(" << (j + 1) << ", "<< prufer[i] << ")  "; 
  				  				 				
                vertex_set[prufer[i] - 1]--; 
  				
				// fill the table of edges
  				edges[l][0]= j + 1;
  				edges[l][1]= prufer[i];
  				l++;
  				
                break; 
            } 
        } 
    } 
  
    j = 0; 
    // For the last element 
    for (int i = 0; i < vertices; i++) 
	{ 
        if (vertex_set[i] == 0 && j == 0) 
		{ 
            //cout << "(" << (i + 1) << ", "; 
            j++;            
            
            // fill the table of edges
  			edges[l][0]= i + 1;
        } 
        else 
			if (vertex_set[i] == 0 && j == 1) 
			{
				//cout << (i + 1) << ")\n"; 
				
				// fill the table of edges
				edges[l][1]= i+1;
			}
    } 
} 

//------------------------------ Sort subtrees ------------------------------
void sort_subtrees(vector<vector<int> >& subtrees, int root,int m)
{
	vector<vector<int> > sorted_subtrees(m);
	int x = 0;
	while((x< subtrees.size())&&(subtrees[x][0]!=root))
	{
		x++;
	}
	sorted_subtrees[0]=subtrees[x];
    subtrees.erase(subtrees.begin()+x);
       	
   	int id=1;
    for(int l=0;l<sorted_subtrees.size();l++)
    {		
       	for(int p=0;p<sorted_subtrees[l].size();p++)
        	for(int q=0;q<subtrees.size();q++)
        		if(subtrees[q][0] == sorted_subtrees[l][p])
        		{
        			sorted_subtrees[id]=subtrees[q];
        			subtrees.erase(subtrees.begin()+q);
        			id++;
        			break;
        		}
	}
	sorted_subtrees.resize(id);
	subtrees.resize(id);
	
	for(int l=id-1;l>=0;l--)
		subtrees[id-1-l]=sorted_subtrees[l];				
}

//-------------- Save distribution items on nodes------------------
void save_dist(int n,int arr[])
{
    char file_name[50];
    char ch[10];
    
	// Get file name
    strcpy(file_name,"nodes_items_");
    sprintf(ch,"%d",n);
    strcat(file_name,ch);
    strcat(file_name,".txt");

    FILE* fichier = NULL;
    fichier = fopen(file_name, "w");
   
    if(fichier!= NULL)
    {
    	for (int i=0;i<n-1;i++)
    		fprintf(fichier, "%d\n", arr[i]);
			
    	fclose(fichier);
    }
    else
    {
        printf("Impossible to open the file %s",file_name);
    }
}


//-------------- Save the tree (set of edges) in a file------------------
void save_tree(int n,int edges[][2])
{
    char file_name[50];
    char ch[10];
    
	// Get file name
    strcpy(file_name,"dependence_graph_");
    sprintf(ch,"%d",n);
    strcat(file_name,ch);
    strcat(file_name,".txt");

    FILE* fichier = NULL;
    fichier = fopen(file_name, "w");
   
    if(fichier!= NULL)
    {
    	for (int i=0;i<n-1;i++)
    		fprintf(fichier, "%d %d %s\n", edges[i][0],edges[i][1],"{}");
			
    	fclose(fichier);
    }
    else
    {
        printf("Impossible to open the file %s",file_name);
    }
}

//-------------------- Create Tree ------------------------
void create_tree(vector<vector<int> >& subtrees,int m,int &root)
{
    // Generate a random fruper sequence
  	int frup[m-2];
    
	// level = 0 that means whatever the level (random)
	int lev;
	//lev = m-2;
	//lev = m/2;
	//lev = m/3;
	lev = m/4;
	fruper(m,lev, frup); 
	
	int i,j;
	int edges[m-1][2]; 
	printTreeEdges(frup, m-2, edges);
	
	// Save the tree into a file
	save_tree(m,edges);
	
	// extract the subtrees and sort them accodring to inclusion in child nodes
	vector<vector<int> > leaves(m);
	int p;
	int k=0,l=-1;
	int ind[m];
	for(i=0;i<m;i++)
	{
		p=0;
		for(j=0;j<m-1;j++)
		{
			if(i+1 == edges[j][1])
			{
				if(p==0)
				{
					l++;
					subtrees[l].push_back(i+1);
					p=1;
				}
				subtrees[l].push_back(edges[j][0]);
			}				
		}
		if(p==0)
		{
			leaves[k].push_back(i+1);
			ind[k]=i+1;
			k++;
		}
	}
	subtrees.resize(l+1);
	
	// Find root
	int vec[m]={0};
	for (int x = 0;  x< subtrees.size(); x++)
        for (int y = 1; y < subtrees[x].size(); y++)
            vec[subtrees[x][y]-1] = 1;
        	
    int x = 0; 
	while(vec[x]!=0)
		x++;
	root = x+1;

	//cout<<"\n\nsubtrees:\n";
	//display_subtrees(subtrees,root);
	
	sort_subtrees(subtrees,root,m);
	
	//cout<<"\nSSSorted subtrees:\n";
	//display_subtrees(subtrees,root);
}
//---------------------------------- Find Leaves of a tree  --------------------------------------- 
void findLeaves(vector<vector<int> > subtrees,int m,vector<int>& leaves)
{
		vector<int> v;
		for(int i=0;i<subtrees.size();i++)
			v.push_back(subtrees[i][0]);
		
		for(int i=1;i<m;i++)
			if(!(find(v.begin(), v.end(),i)!=v.end()))
				leaves.push_back(i);
}
//---------------------------------- Running time --------------------------------------- 
void Display_run_time(vector<vector<vector<double> > > vect_t,int m[])
{  
 	// Save the runinng time for the 4 cases of node number and the 4 cases of knapsack capacity
	vector<vector< vector<double> > > time(4);
	for(int i =0;i<4;i++)
	   	time[i] = vector<vector < double> > (4); // 4 knapsack capacities
	
	for(int i =0;i<4;i++)
    	for(int j =0;j<4;j++)
    			time[i][j] = vector<double>(2); // 2 run time; max and average

	cout<<"\n\n Running times for the 10 instances according to 4 values of Knapsack capcity:\n\n";
	cout<<"Instance|\tk=0\tk=1\tk=2\tk=3\tk=4\tk=5\tk=6\tk=7\tk=8\tk=9\n";
	cout<<"--------------------------------------------------------------------------------------------\n";
	for(int i=0;i<4;i++)
	{
		cout<<"# nodes: "<<m[i];
		cout<<"\n--------\n";
		
		for(int j=0;j<4;j++)
		{
			cout<<"KC ="<<coeff[j]<<"| ";
			for(int k=0;k<10;k++)
				cout<<"\t"<<vect_t[i][j][k];
				cout<<"\n";
			
			time[i][j][0] = max_v(vect_t[i][j]);
			time[i][j][1] = avg_v(vect_t[i][j]);
		}
		cout<<"\n--------------------------------------------------------------------------------------------\n";
	}
	

	cout<<"\n\n";
	cout<<"\t|Sum=0.1\t|Sum=0.2\t|Sum=0.5\t|Sum=0.9\n";
	cout<<"----------------------------------------------------------------------\n";
	cout<<"\t|T_max\tT_avg\t|T_max\tT_avg\t|T_max\tT_avg\t|T_max\tT_avg\n";
	for(int j=0;j<4;j++) // The 4 different vlaues of #nodes
	{
		cout<<"m= "<<m[j]<<"\t";
		for(int l=0;l<4;l++)
		{
			cout<<"|"<<time[j][l][0]<<"\t"<<time[j][l][1]<<"\t";	
		}
		cout<<"\n----------------------------------------------------------------------\n";
	}
	
}

//---------------------------------- Surcharge Operator * for vectors-------------------------------
vector<int> operator*(vector<int> x, vector<int> y)
{
	//cout<<"\n Heeeeere";
	vector<int> r(x.size());
	for (int i=0;i<x.size();i++)
		r[i]=x[i]*y[i];
	return r;
}

// --------------------------- Merging children  --------------------------------------
 vector<int> MergeChildren(vector< vector<int> > WeigChildren, vector< vector<int> > ProfChildren, int W)
{
	vector<int> A(W+1,0); //result vector
	vector<int> H(W+1,0); //needed to save intermediate results
		
	//Initialization of A with the profits of the first child
	for(int i=0;i<WeigChildren[0].size();i++)
		if(WeigChildren[0][i]<=W)
			A[ WeigChildren[0][i] ] = ProfChildren[0][i];
	
	// To ckeck whether at least on item is used  
		
	for(int i=1;i<WeigChildren.size();i++)
	{
		vector<int> index_chang(W+1,0);
		
		for(int j=i;j<=W;j++)
		{
			for(int k=0;k<WeigChildren[i].size();k++)
			{
				int w=j+WeigChildren[i][k]; 
				if (A[j]!=0)
				{
					if(w<=W)
					{
						H[w] = max(H[w],A[j]+ProfChildren[i][k]);
						
						if((A[j]+ProfChildren[i][k])>=H[w])
							index_chang[w]=1;						
					}
				}
				
			}
			H=H*index_chang;
		}
		A=H;
		
	}
	return A;
}


//---------------------------------- DP routine_dependent items ---------------------------------------
int knapSack2_2(int parent,int num_node,vector<vector<int> >& B, int W,int n, vector<int> weig, vector<int> prof)
{

	// Table F for computing function f()	
	vector<vector<int> > F(n+1); 
	for(int i =0;i<n+1;i++)
        F[i] = vector<int> (W+1);
    
    //Initializing F
	for(int i=0;i<=W;i++)
		if(B[parent][i]!=-1)
			F[0][i]=B[parent][i];
		else
			F[0][i] = 0;
    	
    
    // Table new_new for computing function g()	
	vector<int> new_new(W+1); 
	for(int i =0;i<W+1;i++)
        new_new[i] = 0;
	
	for (int i = 1; i <= n; i++) 
	{
    	for (int wt = 1; wt <= W; wt++) 
		{
        	if (weig[i - 1] <= wt)
        	{
        		//cout<<"\n wt = "<<wt<<"\n";
        		if(F[i - 1][wt - weig[i - 1]]>0)
				{
        			F[i][wt] = max(prof[i - 1] + F[i - 1][wt - weig[i - 1]], F[i - 1][wt]);
    				new_new[wt] = max(new_new[wt],prof[i - 1] + F[i - 1][wt - weig[i - 1]]);
				}
				else
        			F[i][wt] = F[i - 1][wt];
			}	        	
        	else
        		F[i][wt] = F[i - 1][wt];
    	}
		
		for(int j =0;j<W+1;j++)
        	B[num_node][j] = new_new[j];
	}
  	
	return F[n][W];
}

//---------------------------------- Deal with the instance i ----------------------------------------------
void deal_instance(int n,int inst,int m,vector<vector<int> > subtrees, int root,int arr[],vector<vector<double> > &vect_t)
{
	// To save execution time (clock)
	double t_run; 
	
	// Get the file name
	char file_name[100]; 
	f_name(file_name,n,inst); 
	
	//fill the list of items for each node from the file
	vector<vector<int> > weights(m);
	vector<vector<int> > profits(m);
	
	int W; // the total sum of the items got from the instance file
	
	load_data_file(file_name, &W, weights, profits,arr);
	
	//Compute the list of knapsack capacity	
	int sum[4];
	list_nbr(sum,coeff,4,W);
	
	
	// Deal with each case of knapsack capacity
	for(int j=1;j<2;j++)  
	{ 
		cout<<"\n\tKnapsack capacity: "<<sum[j];
		vector<vector<int> > B(m, vector<int> (sum[j]+1, 0));
		
		// Initialize B with progits of leaves
		vector<int> leaves;
		findLeaves(subtrees,m,leaves);
		
		for (int i=0;i<leaves.size(); i++)
    	{ 
    		int k = leaves[i];
    		//cout<<"\nheeeere: "<<k;
    		for (int i=0;i<profits[k-1].size();i++)
    			B[k-1][ weights[k-1][i] ] = profits[k-1][i];    					
    	}	
		
		// To save the running time
		cout<<"\n\t\tsubtrees: ";
		
		clock_t tStart = clock();  
		
		// Deal with all subtrees (order of subtrees is inversed)
		for (int k =0;k<subtrees.size(); k++)
    	{	
			cout<<".";
						
			int parent = subtrees[k][0];
			if(subtrees[k].size()>=2)
			{
				if(subtrees[k].size()>2)
				{
					//merging children Cartesian Product		
					vector<int> VW(sum[j]+1);
					for(int w=0;w<=sum[j];w++){VW[w]=w;}
				
					vector<int> merge;
					vector< vector<int> > WeigChildren;
					vector< vector<int> > ProfChildren;
					for(int i=1;i<subtrees[k].size();i++)
					{
						WeigChildren.push_back(VW);		 
						ProfChildren.push_back(B[subtrees[k][i]-1]); 
					}
					
					merge=MergeChildren(WeigChildren, ProfChildren, sum[j]);
					
					B[subtrees[k][1]-1] = merge;				
				}
			
				// Parent is the child node and the child node is the paprent
				int bb = knapSack2_2(subtrees[k][1]-1,parent-1,B,sum[j],weights[parent-1].size(),weights[parent-1],profits[parent-1]);
			}
	
			
    	}
    	
		// Save raunning time    	
    	t_run = (double)(clock() - tStart)/CLOCKS_PER_SEC;
    	vect_t[j][inst]=t_run;
	}
}



//----------------------- Main ----------------------------
int main()
{
	// The number of items: 100,500,1000,5000, or 9000
	int n;
	printf("\n Give the number of items");
	cin>>n; 

    // Save runtime for each 10 execution of the 4 different Knapsack capacities 
	vector<vector<vector<double> > > vect_t(4); 
	for(int i =0;i<4;i++)
        vect_t[i] = vector<vector<double> > (4); 
	for(int i =0;i<4;i++)
		for(int j =0;j<4;j++)
        	vect_t[i][j] = vector<double> (10); 


	// Generate four different # nodes
	int m[4]; 
	list_nbr(m,coeff,4,n);
		
	// Create 4 different trees
	vector< vector<vector<int> > > Trees(4);
	int root[4];
	for(int i=0;i<4;i++)
	{
		Trees[i]=vector<vector<int> > (m[i]);
		create_tree(Trees[i],m[i],root[i]);

		// Display the tree
		display_subtrees(Trees[i],root[i]);
	}
	
	// Deal with the 4 cases of trees
	for(int i=0;i<4;i++)
	{
		cout<<"\n# Nbr nodes = "<<m[i]<<"\n";
		
		//generate randomly the number of items by each node
		int arr[m[i]] = {0}; 
		randomList(m[i], n, arr);
	
		// Deal with the 10 instance
		for(int j=0;j<2;j++)
		{
			cout<<"\n-Instance ="<<j;
			deal_instance(n,j,m[i],Trees[i],root[i],arr,vect_t[i]);
		}
	}
	
	// Display running time				
	Display_run_time(vect_t,m);

	return 0;
}
