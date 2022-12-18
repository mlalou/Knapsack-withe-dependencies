#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <ctime>
#include <cstring>


using namespace std;

int rand_a_b(int a,int b) // generates a random number in [a,b]
{
    int r;
    b; //include b
    if(a>=b)
        r = a;
    else
        r = rand()%(b-a) + a;

    return r;
}

int rand_pos_a_b(int a,int b) // generates a positive random number in [a,b]
{
    int x;
    b; //include b
    if(a>=b)
        x = a;
    else
    {
        do
        {
            x = rand()%(b-a) + a;
        }while(x==0);

    }
    return x;
}



void save_file(int n,int inst)
{
    int i,j,max_gi;
    int r;
    char file_name[50];
    char ch[10];
    double x;

    
	
	// concatenate KP and n
    strcpy(file_name,"KP_");
    sprintf(ch,"%d",n);
    strcat(file_name,ch);
	
	// concatenate BGP_n and number of instance 
    strcat(file_name,"_");
    sprintf(ch,"%d",inst);
    strcat(file_name,ch);
	
	// Add extension "txt"
    strcat(file_name,".txt");

    int w[n]; //contains random values of weights
    int p[n]; //contains random values of profits
    
    for(i=0;i<n;i++)
    {
        w[i]= rand_a_b(1,100);
        

    }
        for(i=0;i<n;i++)
    {
        p[i]= rand_a_b(1,100);
        

    }

    FILE* fichier = NULL;

    fichier = fopen(file_name, "w");
   
    // Calcul de la somme des items
    int s=0; 
    for(i=0;i<n;i++)
    {
        s = s + w[i];
        
    }
    

 	fprintf(fichier, "%d\n", n);
 	fprintf(fichier, "%d\n", s);
 	
 	fprintf(fichier, "\n");

    if(fichier!= NULL)
    {
        for(i=0;i<n;i++)
        {
            
                fprintf(fichier, "%d\n", w[i]);
                fprintf(fichier, "%d\n", p[i]);
                
        }
        fclose(fichier);
    }
    else
    {
        printf("Impossible to open the file %s",file_name);
    }

cout<<"\n saved";

}


//-------------------- Principal --------------------------------
int main()
{
	
	
    int i,nbr;
    int j,k,cc;
    vector< vector<int> > instance;
    int n[5]={100,500,1000,5000,9000};
    srand(time(NULL));

    cout << "\n Give the number of instances for each class : "; 
    cin>>nbr;

    for(i=3;i<5;i++)
    {
    	for(j=0;j<nbr;j++)
            save_file(n[i],j);
            //cout<<"\n File: KP_"<<n[i]<<"_"<<j;
    }

	
    cout<<"\n\n **** FINISH ****";


    return 0;
}
