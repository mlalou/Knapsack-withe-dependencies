#include <bits/stdc++.h>
#include <ctime>
#include<math.h>
#include<stdlib.h>
#include<stdio.h>

#define ll int64_t
#define MX 1000007
#define MX1 1003
using namespace std;

// The instances path
char PATH[100]="E:\\PaperCOR\\1.Experim\\0.Inst\\3.C_stro\\KP_";

//coefficient to compute the number of nodes and sum
float coeff[5]={0.1,0.2,0.5,0.8};

//------------------------- affiche the items filled from the file -----------------------
void display_fillin (vector<vector<int> > weights, vector<vector<int> > profits, int m)
{
	cout<<"\n distribution of items weights and profits:\n";
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<weights[i].size();j++)
			cout<<weights[i][j]<<", "<<profits[i][j]<<"\n";;
		cout<<"\n";	
	}
}
	
//----------------------- print the elements of array, vector and vector of vectors ---------
void printArr(int arr[], int n) 
{ 
    cout<<"\n";
	for (int i = 0; i < n; i++) 
        cout << arr[i] << " "; 
}
void printVect(vector<int> arr) 
{ 
    cout<<"\n";
	for (int i = 0; i < arr.size(); i++) 
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
//----------------------------- print prufer squence ------------------ 
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

//-------------------generate a list of number using coefficients ---------r
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
// --------------------- Compute the max value of a vector -----------------------------------------
double max_v(vector<double> x) 
{
    double m = x[0];
    for(int i=1;i<x.size();i++)
        if(x[i] > m)
            m = x[i];
    return m;
}
//-------------------- Compute the average value of a given vector --------------------
double avg_v(vector<double> x) 
{
    double s,m;
    s=0;
    for(int i=0;i<x.size();i++)
        s = s+ x[i];

    m = s/x.size();
    return m;
}
//---------------------------------- Max of two integer ---------------------------------
int max(int x, int y) {
   return (x > y) ? x : y;
}
//---------------------- Generate m random integers whose sum is n ---------------------- 
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

//----------------------- Get the instance file name ----------------------
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

//--------------------------------- Get the graph name ---------------------
void files_name(char file_name[100],char part_name[100],int n, int inst)
{
	char ch[20];
	
	strcpy(file_name,part_name);
	sprintf(ch,"%d",n);
	strcat(file_name,ch);
	
	// concatenate KP_n and number of instance 
    strcat(file_name,"_");
	sprintf(ch,"%d",inst);
	strcat(file_name,ch);
	
	// Add extension "txt"
    strcat(file_name,".txt");
}

//----------------------------------- Load data from file ---------------------------------------------------
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

//-------------------- Save item distribution among nodes ------------------------
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
	
	cout<<"\nmmmmmmm";
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
		
		cout<<"\nA: "<<endl;
		//for(int i=0;i<A.size();i++)
		//	cout<<" "<<A[i];
	}
	return A;
}

//------------------------------------ Load DAG file ----------------------------
void load_DAG(char file_name[50], vector<vector<int> >& edges)
{	
    FILE* fichier = NULL;
    //vector to save one edge
    vector<int> vect(2);
	char ch[50];
	
	fichier = fopen(file_name, "r");

	if(fichier!= NULL)
     {
     	printf("\n Open successes %s",file_name);
    	do
		{
			fscanf(fichier, "%d %d %s",&vect[0],&vect[1],&ch);		    
			if(vect[0]>=0) //avoid the extra carriage return at the end of the file
				edges.push_back(vect);
			vect[0]=-1; //avoid the extra carriage return at the end of the file
		}      
		while(!feof(fichier));
        fclose(fichier);
    }
    else
    {
        printf("\n Impossible to open the file %s",file_name);
    }	
}
//------------------------------- Load Transitive reduction file ----------------------------
void load_TR(char file_name[100], vector<vector<int> >& TR)
{	
    FILE* fichier = NULL;
    vector<int> vect(2);
	char ch[50];
	
	fichier = fopen(file_name,"r");

	if(fichier!= NULL)
     {
     	printf("\n Open successes %s",file_name);
    	do
		{
			fscanf(fichier, "%d %d %s",&vect[0],&vect[1],&ch);		    
			if(vect[0]>=0) //avoid the extra carriage return at the end of the file
				TR.push_back(vect);
			vect[0]=-1; //avoid the extra carriage return at the end of the file
		}      
		while(!feof(fichier));
        fclose(fichier);
    }
    else
    {
        printf("\n Impossible to open the file %s",file_name);
    }	
}

//----------------------- Load Topological sort file ----------------------------
void load_TS(char file_name[50], vector<int>& TS)
{	
    FILE* fichier = NULL;
	int v;
	
	fichier = fopen(file_name, "r");

	if(fichier!= NULL)
     {
     	printf("\n Open successes %s",file_name);
    	do
		{
			fscanf(fichier, "%d",&v);		    
			TS.push_back(v);
		}      
		while(!feof(fichier));
	
        fclose(fichier);
    }
    else
    {
        printf("\n Impossible to open the file %s",file_name);
    }	
}

//---------------------------- Load Subtrees file ----------------------------
void load_SubT(char file_name[50], vector<vector<int> >& SubT)
{	
    FILE* fichier = NULL;
    vector<int> vect;
	char ch[50];
	
	bool no_digit=true;
	int nbr=0;
	int puis=1;
	fichier = fopen(file_name, "r");
	if(fichier!= NULL)
     {
     	printf("\n Open successes %s",file_name);
    	
    	while(fgets(ch, sizeof(ch), fichier) != NULL) 
		{			
			char *p = ch;
			while (*p) 
			{ 
    			if (isdigit(*p)) 
				{
					no_digit = false;
					long val = strtol(p, &p, 10); 
					nbr=(nbr * puis)+val;
					puis=puis * 10;
    			}	 
				else 
				{
					if (*p == ' ')
					{
						vect.push_back(nbr);
						//printf("\n%d", nbr);						
						nbr = 0;
						puis = 1;
					}
					if (*p == '\n')
					{
						if(no_digit)
						{
						//	vect.push_back(-1);
							SubT.push_back(vect);
							vect.clear();	
						}
						else
						{
							vect.push_back(nbr);
        					SubT.push_back(vect);
							vect.clear();
							//printf("\n%d", nbr);						
							nbr = 0;
							puis = 1;
						}
						no_digit=true;
					}
        			p++;	
    			}
			}
		}
        fclose(fichier);
    }
    else
    {
        printf("\n Impossible to open the file %s",file_name);
    }	
}

//---------------------------- Load Child nodes ----------------------------
void load_ChN(char file_name[50], vector<vector<int> >& ChN)
{	
    FILE* fichier = NULL;
    vector<int> vect;
	char ch[50];
	
	bool no_digit=true;
	int nbr=0;
	int puis=1;
	fichier = fopen(file_name, "r");
	if(fichier!= NULL)
     {
     	printf("\n Open successes %s",file_name);
    	
    	while(fgets(ch, sizeof(ch), fichier) != NULL) 
		{			
			char *p = ch;
			while (*p) 
			{ 
    			if (isdigit(*p)) 
				{
					no_digit = false;
					long val = strtol(p, &p, 10); 
					nbr=(nbr * puis)+val;
					puis=puis * 10;
					//printf("\n%ld", val);						
    			}	 
				else 
				{
					if (*p == ' ')
					{
						vect.push_back(nbr);
						//printf("\n%d", nbr);						
						nbr = 0;
						puis = 1;
					}
					if (*p == '\n')
					{
						if(no_digit)
						{
							//vect.push_back(-1);
							ChN.push_back(vect);
							vect.clear();	
						}
						else
						{
							vect.push_back(nbr);
        					ChN.push_back(vect);
							vect.clear();
							//printf("\n%d",nbr);					
							nbr = 0;
							puis = 1;
						}
						no_digit=true;
					}
        			p++;	
    			}
			}
		}
        fclose(fichier);
    }
    else
    {
        printf("\n Impossible to open the file %s",file_name);
    }	
}

//-------------------------------------- DP routine_dependent items ---------------------------------------
int knapSack2_2(int num_node, int W, vector<vector<int> >& B,vector<int> weig_child,vector<int> prof_child, vector<int> weig,vector<int> prof)
{
	/* weights and profits of the root node
	   are in weig and prof, and those of the 
	   merged children in formal child */
	
	/*
	cout<<"\n";   
	for(int i=0;i<weig.size();i++)
  		cout<<"("<<weig[i]<<","<<prof[i]<<") ";
  		
	cout<<"\ndepend to --->";
		
	for(int i=0;i<weig_child.size();i++)
  		cout<<"("<<weig_child[i]<<","<<prof_child[i]<<") ";
	*/
	
	// Table F for computing function f()	
	vector<vector<int> > F(weig.size()+1); 
	for(int i =0;i<weig.size()+1;i++)
        F[i] = vector<int> (W+1,0);
    
    //Initializing F
	for(int i=0;i<weig_child.size();i++)
		if(weig_child[i]<=W)
			F[0][weig_child[i]]=prof_child[i];
  	
	  cout<<"\n ici ok ok";
  	/*
  	cout<<"\n";
  	for(int i=0;i<=W;i++)
  	cout<<F[0][i]<<" ";
  	*/
    	
    // Table new_new for computing function g()	
	vector<int> new_new(W+1); 
	for(int i =0;i<W+1;i++)
        new_new[i] = 0;
	
	for (int i = 1; i <= weig.size(); i++) 
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
    	
		// printArr(new_new,W+1);
		
		for(int j =0;j<W+1;j++)
        	B[num_node][j] = new_new[j];
	}
  	
  	cout<<"\n ******** F **********";
	/*printMat(F,weig.size(),W);
	
	cout<<"\n\nB["<<num_node+1<<"] = ";
	for(int i=0;i<W+1;i++)
  		cout<<B[num_node][i]<<" ";
	cout<<"\n";	
	
	cout<<"\nF[n][W]= "<<F[weig.size()][W];
	*/
	
	return F[weig.size()-1][W];
}
//----------------------------------------- Predicat for sort() ---------------------------------------
/*bool cmp(const vector<int> &a,const vector<int> &b) 
{ 
	return a.size()>b.size(); 
} 
*/
//---------------------------------- Get intersection between two vectors ----------------------------------
vector<int> intersection(vector<int> &v1, vector<int> &v2)
{
    vector<int> v3;
    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());
    set_intersection(v1.begin(),v1.end(), v2.begin(),v2.end(),back_inserter(v3));
    return v3;
}
//---------------------------------- Get subtrees of the children ----------------------------------
void copy_subtrees(vector<vector<int> > SubT,vector<int> children, vector<vector<int> >& sorted_SubT)
{
	for(int i=0;i<children.size();i++)
		sorted_SubT.push_back(SubT[children[i]]);
}
//---------------------------------- Delete elements of vector from another--------------------------
struct Exist
{
	Exist(const std::vector<int>& vec) : m_vec(vec) {
	}
	bool operator() (int i) {
		return (find(m_vec.begin(), m_vec.end(), i) != m_vec.end());
	}
	const vector<int>& m_vec;
};
//---------------------------------- Check children intersection - form new path ----------------------------------
vector<bool> check_intersection_child(vector<int> children,vector<vector<int> > SubT,vector<int>& formed_path)
{
	
	vector<vector<int> > sorted_SubT;
	copy_subtrees(SubT,children,sorted_SubT);
	
	// We consider only the child node whose subtree has no intresection with the others 
	vector<bool> changed(children.size(),false);
	vector<bool> formed_path_changed;
	
//	sort(sorted_SubT.begin(), sorted_SubT.end(), cmp);
	
//	cout<<"\nsorted subtrees :";
//	printVectVect(sorted_SubT);
		
	for(int i=0;i<children.size()-1;i++)
	{
		for(int j=i+1;j<children.size();j++)
		{	
			// The child to handle is "sorted_SubT[i]" and the rest children are in "sorted_SubT[j]";
			vector<int> intersec;
			
			
			//cout<<"\nChild_SubT: "<<i;
			//printVect(sorted_SubT[i]);
			
			//cout<<"\nRest_Child_SubT: "<<j;
			//printVect(sorted_SubT[j]);
			
			
			// Get intersection between two children
			intersec = intersection(sorted_SubT[i],sorted_SubT[j]);
			
			//cout<<"\nIntersection: ";
			//printVect(intersec);		
			
			if(intersec.size()>0)
			{
				changed[j]=true;		
				// Delete duplicated nodes from other children
				sorted_SubT[j].erase(remove_if(sorted_SubT[j].begin(), sorted_SubT[j].end(), Exist(intersec)), sorted_SubT[j].end());
			}
	
			//cout<<"\nRest_Child_SubT after deletion: ";
			//printVect(sorted_SubT[j]);	
		}
	}
	/*
	cout<<"\n changed : ";
	for(int i=0;i<changed.size();i++)
		cout<<changed[i]<<" ";
	*/
	
		for(int i=0;i<sorted_SubT.size();i++)
		{
			//cout<<"\nSubT size: "<<SubT[children[i]].size();
			if((SubT[children[i]].size()==0)||(changed[i]==true))
			{
				//cout<<"\nch: "<<children[i];
				formed_path.push_back(children[i]); 
				formed_path_changed.push_back(true);
				
				for(int j=0;j<sorted_SubT[i].size();j++)
				{
					formed_path.push_back(sorted_SubT[i][j]);
					formed_path_changed.push_back(true);
				}
			}
			else
			{
				formed_path.push_back(children[i]); 
				formed_path_changed.push_back(false);	
			}
		}
	return formed_path_changed;
}
//---------------------------------- Deal with the instance i ----------------------------------------------
void deal_instance(int n,int inst,int m,vector<vector<int> > subtrees,int arr[],vector<vector<int> > SubT,vector<vector<int> > ChN,vector<int> TS,vector<vector<double> > &vect_t)
{
	// To save execution time (clock)
	double t_run; 
	
	// Get the file name
	char file_name[100]; 
	f_name(file_name,n,inst); 
	
	//fill the list of items for each node from the file
	vector<vector<int> > weights(m);
	vector<vector<int> > profits(m);
	
	int W; // the total sum of the items
	
	/* start */
	load_data_file(file_name, &W, weights, profits,arr); 
	//display_fillin(weights,profits,m);
	//--------
	//char name[50]="test.txt";
	//load_data_file(name, &W, weights, profits,arr);
	/* end */
	
	/* start */	
	//Compute the list of knapsack capacity	
	int sum[4];
	list_nbr(sum,coeff,4,W);
	
	//---------
	//int sum[4];
	//sum[1]=13;
	/* end */
						
	
	/* start */
	// Deal with each case of knapsack capacity
	for(int j=0;j<4;j++)  
	{
		cout<<"\n Knapsack Capacity = "<<sum[j];
	//------------------
	//for(int j=0;j<1;j++)  
	//{
	/* end */

		// To save the running time
		clock_t tStart = clock();  
	
		/* start */
		vector<vector<int> > B(m, vector<int> (sum[j]+1, 0));
		//--------------
		//vector<vector<int> > B(m, vector<int> (sum[1]+1, 0));
		/* end */	
		
		// Deal with all nodes in the topological sort
		for (int k =TS.size()-1;k>=0; k--)
    	{
    		int nd = TS[k];
			cout<<"\n---------------- node: "<<nd<<" --------------";
    		
			// display list of child nodes of TS[k]
    		//cout<<"\nchild nodes:";
			//printVect(ChN[nd]);
    	    	
    		// Check intersection between children and form new path
    		if(ChN[nd].size()>1)
			{
				vector<int> formed_path;
				vector<bool> changed;
				changed = check_intersection_child(ChN[nd],SubT,formed_path);
				
				cout<<"\nFormed path: ";
				printVect(formed_path);	
				
				cout<<"\nchanged : ";
				for(int i=0;i<changed.size();i++)
				cout<<changed[i]<<" ";
			
				/* start */
			
				vector<int> VW(sum[j]+1);
				//-----------
				//vector<int> VW(sum[1]);
				/* end */
			
				
				for(int w=0;w<=sum[j];w++){VW[w]=w;}
				
				
				
				vector<int> merged_child;
				vector< vector<int> > WeigChildren;
				vector< vector<int> > ProfChildren;
				//cout<<"\nbefore merge:";
				for(int i=0;i<formed_path.size();i++)
				{	
					if (changed[i]==true)
					{
						//cout<<"\nhere";
						WeigChildren.push_back(weights[formed_path[i]]); 
						ProfChildren.push_back(profits[formed_path[i]]); 
					}
					else
					{
						WeigChildren.push_back(VW);		 
						ProfChildren.push_back(B[formed_path[i]]); 
					}
					
				}
				
				/*
				cout<<"\n weights: ";
				printVectVect(WeigChildren);
				cout<<"\n pofits: ";
				printVectVect(ProfChildren);
				*/
				
				/* start */
				merged_child=MergeChildren(WeigChildren, ProfChildren, sum[j]);
				//------------------
				//merged_child=MergeChildren(WeigChildren, ProfChildren, sum[1]);
    			/* end */
				
				cout<<"merge done...";
				
				/* start */				
				// Update f()
				int bb = knapSack2_2(nd,sum[j],B,VW,merged_child, weights[nd],profits[nd]);
				//-----------
				//int bb = knapSack2_2(nd,sum[1],B,VW,merged_child, weights[nd],profits[nd]);
				/* end */
    		}
    		else
    		{
    			if((ChN[nd].size()==1)&&(ChN[nd][0]!=-1))
    				/* start */
    				int bb = knapSack2_2(nd,sum[j],B,weights[ChN[nd][0]],profits[ChN[nd][0]], weights[nd],profits[nd]);
    				//------------
    				//int bb = knapSack2_2(nd,sum[1],B,weights[ChN[nd][0]],profits[ChN[nd][0]], weights[nd],profits[nd]);
    				/* end */
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

	// The number of items 100,500,1000,5000, or 9000
	int n;
	printf("\n Give the number of items: ");
	cin>>n; 

    // Save runtime for each 10 execution of the 4 different Knapsack capacities 
	vector<vector<vector<double> > > vect_t(4); 
	for(int i =0;i<4;i++)
        vect_t[i] = vector<vector<double> > (4); 
	for(int i =0;i<4;i++)
		for(int j =0;j<4;j++)
        	vect_t[i][j] = vector<double> (10); 
	
	
	/* start */
	// Generate four different # nodes
	int m[4]; 
	list_nbr(m,coeff,4,n);
	//------------
	//int m[4];
	//m[1]=6;//$$$$$$$$$$$$$$$
	/* end */
	
	/* start */
	// Deal with the 4 cases of DAG
	for(int i=0;i<4;i++)
	{
		cout<<"\n # Nbr nodes = "<<m[i];
	//--------------
	//for(int i=0;i<1;i++)
	//{
	/* end*/
				
		// Load Transitivie reduction of the DAG
		char file_name[100];
		char part_name[100];
		strcpy(part_name,"transitive_reduction_");
		files_name(file_name,part_name,m[i],i);
		/* start */
		//strcpy(file_name,"test.edgelist.txt");
		/* end */
		vector<vector<int> > TR;
		load_TR(file_name,TR);
		
		//printVectVect(TR);
	
	// Load the child nodes of each node
		//char file_name[100];
		strcpy(part_name,"child_nodes_");
		files_name(file_name,part_name,m[i],i);
		/* start */
		//strcpy(file_name,"test.child_nodes.txt");
		/* end */
		vector<vector<int> > ChN;
		load_ChN(file_name,ChN); 
		
		// Load the subtrees of each node
		strcpy(part_name,"subtrees_");
		files_name(file_name,part_name,m[i],i);
		/* start */
		//strcpy(file_name,"test.subtrees.txt");
		/* end */
		vector<vector<int> > SubT;
		load_SubT(file_name,SubT); 
		
		/*
		cout<<"\nChild nodes :";
		printVectVect(ChN);
		cout<<"\nSubtrees :";
		printVectVect(SubT);
		*/
		
		// Load the topological sort
		strcpy(part_name,"topological_sort_");
		files_name(file_name,part_name,m[i],i);
		/* start */
		//strcpy(file_name,"test.topological_sort.txt");//////
		/* end */
		vector<int> TS;
		load_TS(file_name,TS);
	
		
		/* start */
		// generate a random distribution of items among nodes
		int arr[m[i]] = {0}; 
		randomList(m[i], n, arr);
		printArr(arr,m[i]);
		//save_dist(m[i],arr);
		//------------------
		//int arr[6];
		//arr[0]=2;arr[1]=3;arr[2]=1;arr[3]=2;arr[4]=2;arr[5]=2;
		//printArr(arr,6);
		/* end */
		
		/* start */
		// Deal with the 10 instance
		for(int j=0;j<1;j++)
		{
		//-----------------
		//for(int j=0;j<1;j++) 
		//{
		/* end */
			cout<<"\n-------Instance ="<<j;	
			deal_instance(n,j,m[i],TR,arr,SubT,ChN,TS,vect_t[i]);
		}
	}
	
	// Display running time				
	cout<<"\nRunning Times: \n";
	
	Display_run_time(vect_t,m);

	return 0;
}
