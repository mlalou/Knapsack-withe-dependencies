#include <bits/stdc++.h>
#include <ctime>
#include<math.h>
#define ll int64_t
#define MX 1000007
#define MX1 1003
using namespace std;

float coeff[5]={0.1,0.2,0.5,0.8};//coefficient to compute the number of nodes, and the sum

// ---------- Compute the max value in a vector v --------------------
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

    m = s/x.size();
    return m;
}

//----------------------- DP function ------------------- 
bool subset_sum(int root,int parent,int num_node,vector<vector<int> > &B, int sum,int n,vector<int> items)
{

	vector<vector<int> > dp(n+1); 
	for(int i =0;i<n+1;i++)
        dp[i] = vector<int> (sum+1);
	
	//Initializing	
	if(num_node == root)
	{
		for(int i=0;i<=n;i++)
			dp[i][0]=1;
		
		for(int i=1;i<=sum;i++)
    		dp[0][i]=0;
	}
	else
	{
		for(int i=0;i<=n;i++)
			dp[i][0]=0;
		
		for(int i=0;i<=sum;i++)
    		dp[0][i]=B[parent][i];
	}
	//Updating  	
	for(int i=1;i<=n;i++)
	{
    	for(int j=1;j<=sum;j++)
		{
    		if(items[i-1]>j)
        		dp[i][j]=dp[i-1][j];
      		else
      		{
      			dp[i][j]=dp[i-1][j]|dp[i-1][j-items[i-1]];
      			B[num_node][j]= B[num_node][j]|dp[i-1][j-items[i-1]];
			}
        		
    	}
	}
	return dp[n][sum];
}

//----------------------- print the elements of an array ---------
void printArr(int arr[], int n) 
{ 
    cout<<"\n";
	for (int i = 0; i < n; i++) 
        cout << arr[i] << " "; 
} 
  
//---------------------- generate a list of m random non-negative integers whose sum is n 
void randomList(int m, int n,int arr[]) 
{ 
  
    // Create an array of size m where 
    // every element is initialized to 0 
  
    //initialze the m number to 1 to ensure that the random random arr all positive (!=0)
	for (int i = 0; i < m; i++)
    	arr[i]=1;
    	
	srand(time(0)); 
       
	// To make the sum of the final list as n 
    for (int i = 0; i < n-m; i++) { 
  
        // Increment any random element 
        // from the array by 1
        arr[rand() % m]++; 
    } 
  
    // Print the generated list 
    //printArr(arr, m); 
} 

//-------------------- Load data from file -----------------------------
void load_data_file(char file_name[50],int* t_W,vector< vector<int> > &items, int arr[]) 
{
    int i,j;
    FILE* fichier = NULL;
    int x,W;
    
    char ch[20];
	strcpy(ch,"Create_instances/");
    strcat(ch,file_name);
	fichier = fopen(ch, "r");

    if(fichier!= NULL)
    {
        fscanf(fichier, "%d", &x); // the number of nodes
        fscanf(fichier, "%d", &W); // the total weight of all items
        //fscanf(fichier, "%d", &x); // the number of subgraphs
        fscanf(fichier, "%d", &x);
    }

     if(fichier!= NULL)
     {
        i = 0;
        j = 0;
    	do
		{
		    if(j < arr[i])
            {
                items[i].push_back(x); 
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
        printf("Impossible to open the file %s",file_name);
    }
    *t_W=W;
    
}

//------------- Get the file name ----------------
void f_name(char file_name[50],int n, int inst)
{
	char ch[20];
		
	// concatenate KP and n
    strcpy(file_name,"KP_"); 
	sprintf(ch,"%d",n);
	strcat(file_name,ch);
	
	// concatenate KP_n and number of instance 
    strcat(file_name,"_");
	sprintf(ch,"%d",inst);
	strcat(file_name,ch);
	
	// Add extension "txt"
    strcat(file_name,".txt");
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
		// ensure that the first l case have distinct values and then generate random integer (0...l) for the rest of cases
		for (i = 0; i <n-2; i++)
			if(i<l)
				s[i] = i+1;
			else
				s[i] = 1+ (rand() % (l-1));
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

//--------------------- Prints edges of tree represented by a given Prufer code -----------
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
  
    //cout << "\nThe edge set E(G) is :\n"; 
  
    // Find the smallest label not present in 
    // prufer[]. 
    int j = 0; 
    for (int i = 0; i < vertices - 2; i++) { 
        for (j = 0; j < vertices; j++) { 
            // If j+1 is not present in prufer set 
            if (vertex_set[j] == 0) { 
                // Remove from Prufer set and print 
                // pair. 
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

//--------------- This function is used to sort the subtrees of the tree ---------------
bool check_child(vector<vector<int> > ind, int k,vector<int> v2) 
{ 
    // Create a Frequency Table using STL 
    map<int, int> frequency; 
     
    // Increase the frequency of each element 
    // in the frequency table. 
    for (int i = 0; i < k; i++) 
    { 
    	
        frequency[ind[i][0]]++; 
    } 
    
    // Decrease the frequency if the element was found in the frequency table with 
	//the frequency more than 0. else return 0 and if loop is completed return 1. 
    for (int i = 1; i < v2.size(); i++)  
    { 
        if (frequency[v2[i]] > 0) 
            frequency[v2[i]]--; 
        else 
        { 
            return false; 
        } 
    } 
    return true; 
} 

//-------------- display the subtrees of the tree------------------
void display_subtrees(vector<vector<int> > subtrees)
{
	cout<<"Nbr_nodes = "<<subtrees.size()<<"\n";
	for (int i = 0; i < subtrees.size(); i++)
    {
        for (int j = 0; j < subtrees[i].size(); j++)
        {
            cout << subtrees[i][j] << " ";
        }    
        cout << endl;
    }
}

//-------------------- Create Tree ------------------------

void create_tree(vector<vector<int> >& new_sorted_subtrees,int m,int &root)
{
  	
    // Generate a random fruper sequence
  	int frup[m-2];

    //cout<<"\n level ="<<level;
	fruper(m, 0, frup); // level = 0 that means whatever the level (random)
	
	//print_fruper(frup,m-2);
	
	int i,j;
	int edges[m-1][2]; 
	printTreeEdges(frup, m-2, edges);
	
    // extract the subtrees and sort them accodring to inclusion in child nodes
	vector<vector<int> > subtrees(m);
	vector<vector<int> > sorted_subtrees(m);
	vector<vector<int> > leaves(m);
	int p;
	int k=0,l=-1;
	int ind[m];
	for(i=0;i<m;i++)
	{
		//cout<<"\n ---- m = "<<i;
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
	 	
	display_subtrees(subtrees);
	
	//------------------- Find root ---------------- 

	int vec[m]={0};
	for (int x = 0;  x< subtrees.size(); x++)
        for (int y = 1; y < subtrees[x].size(); y++)
            vec[subtrees[x][y]-1] = 1;
        	
    int x = 0; 
	while(vec[x]!=0)
		x++;
	root = x+1;
	cout<<"\n ---------- root ="<<root;
	//--------------------------------------------------
	
	x = 0;
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
	new_sorted_subtrees.resize(id);
	
	cout<<"\n\nsorted subtrees:\n";
	display_subtrees(sorted_subtrees);
	
	for(int l=id-1;l>=0;l--)
		new_sorted_subtrees[id-1-l]=sorted_subtrees[l];				
	
	display_subtrees(new_sorted_subtrees);
		
		
  
	
}

//-------------------- Deal with the instance i -----------
void deal_instance(int n,int inst,int m,vector<vector<int> > subtrees, int root,vector<vector<double> > &vect_t)
{
	double t_run; //save execution time (clock)
	
	//get the file name
	char file_name[50]; 
	f_name(file_name,n,inst); 
	//f_name(file_name,n,0); // ici 
	cout<<"\nFile name: "<<file_name;
		
	//generate randomly the number of items by each node
	int arr[m] = {0}; 
	randomList(m, n, arr);	
		
	//fill the list of items for each node from the file
	vector< vector<int> > items(m);
	int W;
	
	load_data_file(file_name, &W, items,arr);
					
	//Compute the list of sum	
	int sum[4];
	list_nbr(sum,coeff,5,W);
	
	//deal with each case of sum
	for(int j=0;j<4;j++)
	{
		cout<<"\nSum = "<<sum[j];
		
		vector<vector<int> > B(m, vector<int> (sum[j]+1, 0));
		
		int cpt=0;//compteur juste pour compter le nomber de sous-arbre traité
		
		clock_t tStart = clock();  // Used to save the running time
		
		for (int k = subtrees.size()-1; k <=subtrees.size()-1 ; k--)// deal with all subtrees (because the order of subtrees is inverse (>=0))
    	{	
			int parent = subtrees[k][0];
					    		
			for (int l = 0; l <subtrees[k].size() ; l++) // deal with all nodes of the considered subtree
       		{
       			
       			if((k != subtrees.size()-1)&&(l==0)) // if it is not the first subtree for which the root is the root of the whole tree, then the root is already treated and so skip it
				   		continue;
				   		
       			int bb = subset_sum(root-1,parent-1,subtrees[k][l]-1,B,sum[j],items[subtrees[k][l]-1].size(),items[subtrees[k][l]-1]);
			}           		
    	}
    	
    	t_run = (double)(clock() - tStart)/CLOCKS_PER_SEC;
    	vect_t[j][inst]=t_run;
	}
	
}

//----------------------- Main ----------------------------
int main()
{
	int n[5]={100,500,1000,5000,9000}; //the number of items
	int m[4]; // the number of nodes
	
	vector< vector<vector< vector<double> > > > time(5); // 5 classes d'instances
	
	for(int i =0;i<5;i++)
	   	time[i] = vector< vector<vector < double> > >(4); // 4 numbers of nodes
	
	for(int i =0;i<5;i++)
    	for(int j =0;j<4;j++)
    			time[i][j] = vector<vector < double> > (4); // 4 valeurs de sommes différentes

    for(int i =0;i<5;i++)
    	for(int j =0;j<4;j++)
    		for(int k =0;k<4;k++)
        			time[i][j][k] = vector<double>(2); // il y a 2 temps max et average
	
	
	//créer le vecteur pour sauvegarder les temps d'exécution des 10 instance, on a pour chaque sum Tmax et Tavg (voir paper) qu'on calcule des 10 exécutions, donc 10 temps pour chaque sum
	vector<vector<double> > vect_t(4);
	for(int i =0;i<4;i++)
        vect_t[i] = vector<double> (10); 

	for(int i=4;i<5;i++)//deal with all classes of instances (<5)
	{
		cout<<"\nInstance - "<<n[i]<<" items:";
		
		//generate the number of nodes
		list_nbr(m,coeff,4,n[i]);
		
		//m[0]=10;
		
		for(int j=0;j<4;j++)//deal with all case of number of nodes m (<4)
		{
			cout<<"\nNumber_nodes : "<<m[j];
			
			//generate the tree
			int root;
			vector<vector<int> > subtrees(m[j]);
			create_tree(subtrees,m[j],root);
				
			for(int k=0;k<10;k++)//deal with the 10 instance of each case (<10)
				deal_instance(n[i],k,m[j],subtrees,root,vect_t);

			// ---------------------affichage des temps d'exécution -----------------------------------------------------
			
			cout<<"\n\nRunning times for the 10 instances of a given class according to 4 values of Sum:\n\n";
			cout<<"Instance|\tk=0\tk=1\tk=2\tk=3\tk=4\tk=5\tk=6\tk=7\tk=8\tk=9\n";
			cout<<"--------------------------------------------------------------------------------------------\n";
			for(int p=0;p<vect_t.size();p++)
			{
				cout<<"Sum ="<<coeff[p]<<"| ";
				for(int j=0;j<vect_t[p].size();j++)
					cout<<"\t"<<vect_t[p][j];
				cout<<"\n";
			}
			for(int l=0;l<4;l++)
			{
				time[i][j][l][0] = max_v(vect_t[l]);
				time[i][j][l][1] = avg_v(vect_t[l]);
			}
		}
		
		cout<<"\nInstance_Classe "<<n[i]<<"\n";
		
		cout<<"\t|Sum=0.1\t|Sum=0.2\t|Sum=0.5\t|Sum=0.9\n";
		cout<<"----------------------------------------------------------------------\n";
		cout<<"\t|T_max\tT_avg\t|T_max\tT_avg\t|T_max\tT_avg\t|T_max\tT_avg\n";
		for(int j=0;j<4;j++) //all case of number of nodes m (<4)
		{
			cout<<"m= "<<m[j]<<"\t";
			for(int l=0;l<4;l++)
			{
				cout<<"|"<<time[i][j][l][0]<<"\t"<<time[i][j][l][1]<<"\t";	
			}
			cout<<"\n----------------------------------------------------------------------\n";
		}	
	}
  	return 0;
}
