#include <iostream.h>
#include <fstream.h>
#include <math.h>
typedef char boolean;
const int numSaveMechs = 4; // number of mechanisms required
const int xRestr = 1, yRestr = 2, rotnRestr = 4; // restraint codes
const int RESTRAINED = -1; // used to mark restrained degrees of freedom
const int maxNodes = 60, maxElements = 60, maxNumDOF = 180; // problem size constants
const double VERYLARGE = 1E12, PRECISION = 1E-12; // other constants

struct nodeRecord {
  double x,y,load[2]; // first load is in x dirn, second in y
  int restr,dof[2]; // restraint code, DOF numbers
  int numConnections, connection[6]; // connecting members (max of 6)
};

struct elementRecord {
  int node[2];
  double plasticMoment;
};

class mechanismClass {
  public:
    double *disp, *rotn; // joint displacements and member rotations
    boolean done; // working flag
    
    mechanismClass(); // constructor
    ~mechanismClass() {delete disp; delete rotn;}; // destructor
    void calcRotations(); // calculates rotations from displacements
    double calcRelativeRotn(int m1, int m2) {
      return rotn[m1]-rotn[m2];
    };
    double calcCollapseLoad() {return calcIntVW()/calcExtVW();};
  private:
    double calcExtVW();
    double calcIntVW();
};

class combinedMechClass : public mechanismClass {
  public:
    double *factor; // factors for independent mechanisms
    boolean *mechIncluded; // flags for mechanisms included
    combinedMechClass(); // constructor
    ~combinedMechClass(); // destructor
    void combine(double lastLoad);
    void addDisps(int n);
    void addMech(int i, double rotn1, double rotn2, double lastLoad);
  private: // used to ensure same factor isn't repeated plastic
    int numUsed;
    double *usedFactor;
    void resetDone() {numUsed = 0;};
    boolean isDone(double f);
    void setDone(double f) {usedFactor[numUsed++] = f;};
};

/* function prototypes */
void ReadData();
void AssignDOF();
void FormIndependentMechanisms();
void ComputeCollapseMechanism();
void PrintResults();
void SortCollapseLoads();

boolean NotInList(double factor, const combinedMechClass& mech);

void Swap(double& a, double& b) {double temp = a; a = b; b = temp;}
/* global data */
// specified data
int numNodes,numElements;
elementRecord element[maxElements];
nodeRecord node[maxNodes];

// working variables
int numDOF,numIndepMechs,maxHinges;
mechanismClass *independentMech;
long mechCount = 0L;
// results
double collapseLoadFactor[numSaveMechs];
combinedMechClass *collapseMechanism[numSaveMechs];
/* main program */
void main() {
  ReadData();
  AssignDOF();
  FormIndependentMechanisms();
  ComputeCollapseMechanism();
  PrintResults();
}
/* independent functions */
void ReadData() {
  char fileName[16];
  cout << "Enter file name: ";
  cin >> fileName;
  ifstream infile(fileName);
  infile >> numNodes;
  for(int i=0;i<numNodes;i++)
  infile>>node[i].x>>node[i].y>>node[i].restr>>node[i].load[0]>>node[i].load[1];
  infile >> numElements;
  for(int i=0;i<numElements;i++)
  infile>>element[i].node[0]>>element[i].node[1]>>element[i].plasticMoment;
} // end of ReadData()

void AssignDOF() {
  numDOF = maxHinges = 0;
  for(int i=0;i<numNodes;i++) { // Alle Knoten
    node[i].dof[0] = node[i].dof[1] = RESTRAINED;
    if(!(node[i].restr&xRestr)) node[i].dof[0] = numDOF++;
    if(!(node[i].restr&yRestr)) node[i].dof[1] = numDOF++;
    node[i].numConnections = 0;
  }
  for(int i=0;i<numElements;i++)
    for(int j=0;j<2;j++) {
    int n = element[i].node[j];
    node[n].connection[node[n].numConnections++] = i;
    maxHinges++;
  }
} // end of AssignDOF()
void FormIndependentMechanisms() {
/* -------------------- Form compatibility matrix ----------------------------- */
  double C[maxElements][maxNumDOF]; // should be dynamically allocated
  for(int i=0;i<numElements;i++) {
    for(int j=0;j<numDOF;j++) C[i][j] = 0.0;
      double a[4]; // compatibility of element stretching with global DOF
      double dx = node[element[i].node[1]].x - node[element[i].node[0]].x;
      double dy = node[element[i].node[1]].y - node[element[i].node[0]].y;
      double len = sqrt(dx*dx + dy*dy);
      a[0] = -dx/len, a[1] = -dy/len, a[2] = dx/len, a[3] = dy/len;
  /* assemble element compatibility into global compatibility */
    for(int j0=0;j0<2;j0++)
      for(int k0=0;k0<2;k0++) {
        int n0 = node[element[i].node[j0]].dof[k0];
        if(n0!=RESTRAINED) C[i][n0] = a[j0*2 + k0];
      }
    }
  /* -------------------- Perform Gaussian elimination -------------------------- */
  numIndepMechs = 0;
  int *dependentCol = new int [numDOF]; // numIndepMechs can't be more than numDOF
  int k = 0;
  for(int i=0;i<numElements;i++) {
    while(k<numDOF && C[i][k]==0.0) {
      int j = i+1;
      while(j<numElements && C[j][k]==0.0) j++;
        if(j!=numElements) // swap rows i and j
          for(int j1=0;j1<numDOF;j1++)
            Swap(C[j][j1],C[i][j1]);
        else dependentCol[numIndepMechs++] = k++;
    }
    double d = C[i][k];
    for(int j=k;j<numDOF;j++) C[i][j] /= d;
      for(int i1=0;i1<numElements;i1++)
        if(i1!=i) {
          double d=C[i1][k];
          for(int j=k;j<numDOF;j++) C[i1][j] -= C[i][j]*d;
        }
    k++;
  }
  while(k<numDOF) dependentCol[numIndepMechs++] = k++;
  /* -------------------- Extract independent mechanisms -------------------------- */
  ::independentMech = new mechanismClass [numIndepMechs];
  k=0;
  int j=0;
  for(int i=0;i<numElements;i++) {
    while(dependentCol[j]==k) independentMech[j++].disp[k++] = -1.0;
    for(int i1=j;i1<numIndepMechs;i1++)
      independentMech[i1].disp[k] = C[i][dependentCol[i1]];
    k++;
  }
  for(int i=0;i<numIndepMechs;i++) {
    independentMech[i].calcRotations();
    independentMech[i].done = false;
  }
  delete [] dependentCol; // clean up
  /* --------- Now 'purify' mechanisms so that none is contained by another -------- */
  boolean allDone = false;
  while(!allDone) {
    allDone = true;
    for(int i=0;i<numIndepMechs;i++)
      if(!independentMech[i].done) {
        for(int j=0;j<numIndepMechs;j++)
          if(j!=i) {
            boolean contained = true, first = true;
            double ratio;
            for(int k=0;k<numNodes;k++) {
              // only consider nodes which move in this mechanism
              int m = node[k].dof[0];
              if(m!=RESTRAINED && fabs(independentMech[i].disp[m])<PRECISION) {
                m = node[k].dof[1];
                if(m!=RESTRAINED && fabs(independentMech[i].disp[m])<PRECISION)
                  continue; // skip this node
              }
              // check each possible hinge location at this node
              for(int mi=0;mi<node[k].numConnections;mi++) {
                int m1 = node[k].connection[mi];
                double rotni = independentMech[i].rotn[m1];
                double rotnj = independentMech[j].rotn[m1];
                for(int mj=mi+1;mj<node[k].numConnections;mj++) {
                  int m2 = node[k].connection[mj];
                  double relRotni = rotni - independentMech[i].rotn[m2];
                  if(fabs(relRotni)>PRECISION) {
                    double relRotnj = rotnj - independentMech[j].rotn[m2];
                    if(fabs(relRotnj)<PRECISION) {
                      contained = false; break;
                    } else if(first) ratio = relRotnj/relRotni, first = false;
                    else if(fabs(ratio-relRotnj/relRotni)>PRECISION){
                      contained = false; 
                      break;
                    }
                  }
                }
              }
            }
            // remove mechanism i from mechanism j if contained
            if(!first && contained) {
              for(k=0;k<numElements;k++)
                independentMech[j].rotn[k] -= ratio*independentMech[i].rotn[k];
              for(k=0;k<numDOF;k++)
                independentMech[j].disp[k] -= ratio*independentMech[i].disp[k];
              independentMech[j].done = false;
              allDone = false;
            }
          }
        independentMech[i].done = true;
      }
  }
} // end of FormIndependentMechanisms()

void ComputeCollapseMechanism() {
for(int i=0;i<numSaveMechs;i++) {
collapseMechanism[i] = new combinedMechClass;
collapseLoadFactor[i] = VERYLARGE;
}
for(int i=0;i<numIndepMechs;i++) {
combinedMechClass tempMech;
tempMech.factor[i] = 1.0;
tempMech.mechIncluded[i] = true;
tempMech.addDisps(i);
tempMech.combine(VERYLARGE);
}
} // end of ComputeCollapseMechanism()
void PrintResults() {
for(int k=0;k<numSaveMechs;k++) {
cout << "\n\nCollapse load factor " << (k+1) << " is " <<
collapseLoadFactor[k];
cout << "\n member rotations: ";
for(int i=0;i<numElements;i++) {
double rotn = 0.0;
for(int j=0;j<numIndepMechs;j++)
rotn +=
collapseMechanism[k]->factor[j]*independentMech[j].rotn[i];
cout << rotn << " ";
}
cout << "\n joint displacements: ";
for(int i=0;i<numNodes;i++)
for(int j=0;j<2;j++) {
double disp = 0.0;
int n = node[i].dof[j];
if(n!=RESTRAINED)
for(int k1=0;k1<numIndepMechs;k1++)
disp +=
collapseMechanism[k]->factor[k1]*independentMech[k1].disp[n];
cout << disp << " ";
}
}
cout << "\n\nNumber of mechanisms processed was " << mechCount << "\n";
} // end of PrintResults()
void SortCollapseLoads() { // bubble sort - only a few to sort
boolean changed;
do{
changed = false;
for(int j=0;j<numSaveMechs-1;j++)
if(collapseLoadFactor[j+1]<collapseLoadFactor[j]) {
changed = true;
Swap(collapseLoadFactor[j],collapseLoadFactor[j+1]);
for(int k=0;k<numIndepMechs;k++)
Swap(collapseMechanism[j]->factor[k],collapseMechanism[j+1]->factor[k]);
}
} while(changed);
} // end of SortCollapseLoads()
boolean NotInList(double factor, const combinedMechClass& mech) {
for(int i=0;i<numSaveMechs;i++)
if(fabs(factor-collapseLoadFactor[i])<PRECISION) {
boolean same = true, first = true;
double ratio = 0.0;
for(int j=0;j<numIndepMechs;j++)
if(fabs(mech.factor[j])>PRECISION)
if(first)
ratio =
collapseMechanism[i]->factor[j]/mech.factor[j], first = false;
else
if(fabs(ratio-collapseMechanism[i]->factor[j]/mech.factor[j])>PRECISION)
{same = false; break;}
if(same) return false;
}
return true;
} // end of NotInList(...)
/*---------------------- mechanismClass implementation
--------------------------*/
mechanismClass::mechanismClass() {
disp = new double [numDOF]; for(int i=0;i<numDOF;i++) disp[i] =
0.0;
rotn = new double [numElements]; for(int i=0;i<numElements;i++) rotn[i]
= 0.0;
}
void mechanismClass::calcRotations() {
for(int i=0;i<numElements;i++) {
double dx = node[element[i].node[1]].x -
node[element[i].node[0]].x;
double dy = node[element[i].node[1]].y -
node[element[i].node[0]].y;
double len2 = dx*dx + dy*dy;
double dd[2];
for(int j=0;j<2;j++) {
dd[j] = 0.0;
int r0 = node[element[i].node[0]].dof[j];
if(r0!=RESTRAINED) dd[j] -= disp[r0];
int r1 = node[element[i].node[1]].dof[j];
if(r1!=RESTRAINED) dd[j] += disp[r1];
}
rotn[i] = (-dd[0]*dy + dd[1]*dx)/len2;
}
}
double mechanismClass::calcExtVW() {
double extVW = 0.0;
for(int i=0;i<numNodes;i++)
for(int j=0;j<2;j++) {
int r = node[i].dof[j];
if(r>=0) extVW += node[i].load[j]*disp[r];
}
return extVW;
}
double mechanismClass::calcIntVW() {
double intVW = 0.0;
for(int i=0;i<numNodes;i++)
if(node[i].restr&rotnRestr) // joint has zero rotn, hinges form
in all members
for(int j=0;j<node[i].numConnections;j++) {
int m = node[i].connection[j];
intVW +=
fabs(rotn[m]*element[m].plasticMoment);
}
else { // must find out which hinge mechanism leads to lowest
internal work
double jointIntVW = VERYLARGE;
for(int j=0;j<node[i].numConnections;j++) {
double jointRotn = rotn[node[i].connection[j]];
double thisVW = 0.0;
for(int k=0;k<node[i].numConnections;k++) {
int m = node[i].connection[k];
thisVW += fabs((rotn[m] -
jointRotn)*element[m].plasticMoment);
}
if(thisVW<jointIntVW) jointIntVW = thisVW;
}
intVW += jointIntVW;
}
return intVW;
}
/*---------------------- combinedMechClass implementation
------------------------*/
combinedMechClass::combinedMechClass() {
factor = new double [::numIndepMechs];
mechIncluded = new boolean [::numIndepMechs];
for(int i=0;i<numIndepMechs;i++) {
factor[i] = 0.0;
mechIncluded[i] = false;
}
usedFactor = new double [maxHinges];
}
combinedMechClass::~combinedMechClass()
{delete factor; delete usedFactor; delete mechIncluded;}
void combinedMechClass::combine(double lastLoad) {
// check the collapse load of this mechanism
double collapse = calcCollapseLoad();
double absCollapse = fabs(collapse);
if(absCollapse>lastLoad && lastLoad!=VERYLARGE) return;
mechCount++;
if(absCollapse<collapseLoadFactor[numSaveMechs-1]) { // save this one
if(collapse<0.0) { // change sign
collapse = -collapse;
for(int i=0;i<numIndepMechs;i++) factor[i] =
-factor[i];
for(int i=0;i<numDOF;i++) disp[i] = -disp[i];
for(int i=0;i<numElements;i++) rotn[i] = -rotn[i];
}
if(NotInList(collapse,*this)) { // haven't got this one yet
int n = numSaveMechs - 1;
collapseLoadFactor[n] = collapse;
for(int i=0;i<numIndepMechs;i++)
collapseMechanism[n]->factor[i] = factor[i];
SortCollapseLoads();
}
}
// now find all mechanisms which can be formed by the addition of one
more
for(int i=0;i<numIndepMechs;i++)
if(!mechIncluded[i]) {
resetDone();
for(int nn=0;nn<numNodes;nn++)
for(int
nci=0;nci<node[nn].numConnections;nci++) {
int m1 = node[nn].connection[nci];
if(node[nn].restr&rotnRestr) { // joint
has zero rotation
double rotn1 = rotn[m1];
double rotn2 =
independentMech[i].rotn[m1];
addMech(i,rotn1,rotn2,absCollapse);
}
else
for(int
ncj=nci+1;ncj<node[nn].numConnections;ncj++) {
int m2 =
node[nn].connection[ncj];
double rotn1 =
calcRelativeRotn(m1,m2);
double rotn2 =
independentMech[i].calcRelativeRotn(m1,m2);
addMech(i,rotn1,rotn2,absCollapse);
}
}
}
}
void combinedMechClass::addMech(int n, double rotn1, double rotn2, double
lastLoad) {
if(fabs(rotn1)>PRECISION && fabs(rotn2)>PRECISION) {
double f = -rotn1/rotn2;
if(!isDone(f)) {
mechIncluded[n] = true;
combinedMechClass tempMech;
for(int j=0;j<numIndepMechs;j++)
{tempMech.factor[j] = factor[j];
tempMech.mechIncluded[j] = mechIncluded[j];}
for(int j=0;j<numDOF;j++) tempMech.disp[j] = disp[j];
for(int j=0;j<numElements;j++) tempMech.rotn[j] =
rotn[j];
tempMech.factor[n] = -rotn1/rotn2;
tempMech.addDisps(n);
tempMech.combine(lastLoad);
setDone(f);
}
}
}
void combinedMechClass::addDisps(int n) {
if(factor[n]!=0.0) {
for(int i=0;i<numDOF;i++) disp[i] +=
factor[n]*independentMech[n].disp[i];
for(int i=0;i<numElements;i++) rotn[i] +=
factor[n]*independentMech[n].rotn[i];
}
}
boolean combinedMechClass::isDone(double f) {
for(int i=0;i<numUsed;i++)
if(fabs(usedFactor[i]-f)<PRECISION)
return true;
return false;
}
/*--------------------------------- end of program
------------------------------------*/