#include <iostream>
//#include <fstream>
//#include <sstream>
#include <cmath>
#include "apmatrix.h"
#include <time.h>
#include <chrono>
#include "windows.h"

void printUM(apmatrix<double> unmod);
void printM(apmatrix<double> matrix);

apmatrix<double>add(apmatrix<double>mtx, apmatrix<double>matrix);
apmatrix<double>subtract(apmatrix<double>mtx, apmatrix<double> matrix);
apmatrix<double>mult(apmatrix<double>mtx, apmatrix<double> matrix);
apmatrix<double>rE(apmatrix<double> matrix, apmatrix<double>unmod);
apmatrix<double>rrE(apmatrix<double> matrix, apmatrix<double>unmod);
apmatrix<double>transpose(apmatrix<double>matrix, apmatrix<double>unmod);
apmatrix<double>inverse(apmatrix<double> mtx, apmatrix<double>unmod);
apmatrix<double>linReg(apmatrix<double>x,apmatrix<double>y);
apmatrix<double>regrssn(apmatrix<double>x,apmatrix<double>y, apmatrix<string> f, apmatrix<double> v);
apmatrix<double>augment(apmatrix<double>A,apmatrix<double>B);
apmatrix<double>redaug(apmatrix<double>aug, apmatrix<double> A, apmatrix<double>B);
apmatrix<double>proj(apmatrix<double>v, apmatrix<double>u);
apmatrix<double>orthogonal(apmatrix<double>A);
apmatrix<double>orthonormal(apmatrix<double>V);
apmatrix<double> eV(apmatrix<double>A, apmatrix<double> Q, apmatrix<double>R1);
void enter();
void setMtx(int indx, apmatrix<double>mtx);
void store(apmatrix<double>matrix,apmatrix<double>orig);
void rVals();
//void onOpen();
//void onClose();
apmatrix<double> getMtx(int val);
apmatrix<double> getUnM(int val);
apmatrix<double> cmatrix,m1,m2,m3,m4,m5,m6,unmod,unmod1, unmod2,unmod3,unmod4,unmod5,unmod6, dataX, dataY, dataV;
apmatrix<string> dataF;
int row, col;
double spent;
using namespace std;

int main(int argc, char *argv[])
{
    SetConsoleTitle("Matrix Operations");
    bool menu=true, loop = false, tp = false, red = false;
//cout<<"Opened";
// onOpen();
    while(menu) //Menu
    {
        cout<<"\nMenu\n"<<"----\n";
        cout<<"1) Enter Matrix\n"<<"2) Show Matrix\n"<<"3) Transpose Matrix\n"<<"4) Add Matrices\n"<<"5) Subtract Matrices\n"<<"6) Multiply Matrices\n"<<"7) Row Echelon\n"<<"8) Reduced Row Echelon\n"<<"9) Inverse Matrix\n"<<"10) Linear Regression\n"<<"11) Matrix Regression\n"<<"12) QR Factorization\n"<<"0) Quit\n"<<"----\n>";
        int val;
        cin>>val;
        if(val == 1)
        {
            enter();
            system("PAUSE");
        }
        else if(val==2)   //Show Matrix
        {
            bool l;
            while(l)
            {
                cout<<"Which matrix do you want to see (1-6)? "<<"\n>";
                int show;
                cin>>show;
                if(show>=1 && show<=6)
                {
                    printM(getMtx(show));
                    l=false;
                }
                else
                {
                    cout<<"Please enter a positive integer from 1-6.\n";
                }
            }
            cout<<"\n";
            system("PAUSE");
        }
        else if(val == 3) //Transpose Matrix
        {
            int num;
            cout<<"Which matrix do you want to transpose?"<<"\n>";
            cin>>num;
            //Changes original matrix to matrix stored in num value.
            setMtx(num,getMtx(num));
            cmatrix = transpose(getMtx(num), getUnM(num));
            printUM(getUnM(num));
            cout<<"Transposed Matrix:"<<"\n";
            printM(cmatrix);
            store(cmatrix, getUnM(num));
            system("PAUSE");
        }
        else if(val == 4) //Add Matrices
        {
            int mtx1, mtx2;
            apmatrix<double>sum;
            cout<<"Enter first matrix to add: ";
            cin>>mtx1;
            cout<<"Enter second matrix to add: ";
            cin>>mtx2;
            sum = add(getMtx(mtx1), getMtx(mtx2));
            if(sum[1][1] == 1)
            {
                cout<<"Matrix Size Mismatch. <"<<getMtx(mtx1).numrows()<<","<<getMtx(mtx1).numcols()<<"> + <"<<getMtx(mtx2).numrows()<<","<<getMtx(mtx2).numcols()<<"> is not possible.\n";
                system("PAUSE");
            }
            else
            {
                printM(getMtx(mtx1));
                //cout<<"\n";
                for(int i=0; i<getMtx(mtx1).numcols()/2; i++)
                    cout<<" ";
                cout<<"+\n";
                printM(getMtx(mtx2));
                cout<<"\n-";
                for(int i=0; i<sum.numcols(); i++)
                    cout<<"-";
                cout<<"\n";
                printM(sum);
                store(sum, sum);
            }
            system("PAUSE");
        }
        else if(val == 5) //Subtract Matrices
        {
            int mtx1, mtx2;
            apmatrix<double>diff;
            cout<<"Enter first matrix to subtract: ";
            cin>>mtx1;
            cout<<"Enter second matrix to subtract: ";
            cin>>mtx2;
            diff = subtract(getMtx(mtx1), getMtx(mtx2));
            if(diff[1][1] == 1)
            {
                cout<<"Matrix Size Mismatch. <"<<getMtx(mtx1).numrows()<<","<<getMtx(mtx1).numcols()<<"> - <"<<getMtx(mtx2).numrows()<<","<<getMtx(mtx2).numcols()<<"> is not possible.\n";
                system("PAUSE");
            }
            else
            {
                printM(getMtx(mtx1));
                //cout<<"\n";
                for(int i=0; i<getMtx(mtx1).numcols()/2; i++)
                    cout<<" ";
                cout<<"-\n";
                printM(getMtx(mtx2));
                cout<<"\n-";
                for(int i=0; i<diff.numcols(); i++)
                    cout<<"-";
                cout<<"\n";
                printM(diff);
                store(diff, diff);
            }
            system("PAUSE");
        }
        else if(val == 6) //Multiply Matrices
        {
            int mtx1, mtx2;
            apmatrix<double>prod;
            cout<<"Enter first matrix to multiply: ";
            cin>>mtx1;
            cout<<"Enter second matrix to multiply: ";
            cin>>mtx2;
            prod = mult(getMtx(mtx1), getMtx(mtx2));
            if(prod[1][1] == 1)
            {
                cout<<"Matrix Size Mismatch. <"<<getMtx(mtx1).numrows()<<","<<getMtx(mtx1).numcols()<<"> x <"<<getMtx(mtx2).numrows()<<","<<getMtx(mtx2).numcols()<<"> is not possible.\n";
                system("PAUSE");
            }
            else
            {
                printM(getMtx(mtx1));
                //cout<<"\n";
                for(int i=0; i<getMtx(mtx1).numcols()/2; i++)
                    cout<<" ";
                cout<<"*\n";
                printM(getMtx(mtx2));
                cout<<"\n-";
                for(int i=0; i<prod.numcols(); i++)
                    cout<<"-";
                cout<<"\n";
                printM(prod);
                store(prod, prod);
            }
            system("PAUSE");
        }
        else if(val == 7) //Row Echelon
        {
            int num;
            cout<<"Which matrix do you want to reduce?\n>";
            cin>>num;
            //Changes original/unmod matrix to matrix stored in num value.
            setMtx(num,getMtx(num));
            printUM(getUnM(num));
            cmatrix.resize(getMtx(num).numrows(), getMtx(num).numcols()); //resize to fit target matrix
            cmatrix = rE(getMtx(num), getUnM(num)); //matrix = target matrix
            cout<<"REF Matrix:"<<"\n";
            printM(cmatrix); //Print new matrix
            store(cmatrix,getUnM(num)); //Store matrix
            system("PAUSE");
        }
        else if(val == 8) //Reduced Row Echelon
        {
            int num;
            cout<<"Which matrix do you want to reduce?\n>";
            cin>>num;
            //Changes original/unmod matrix to matrix stored in num value.
            setMtx(num,getMtx(num));
            printUM(getUnM(num));
            cmatrix.resize(getMtx(num).numrows(), getMtx(num).numcols()); //resize to fit target matrix
            cmatrix = rrE(rE(getMtx(num), getUnM(num)),rE(getMtx(num), getUnM(num))); //matrix = target matrix
            cout<<"RREF Matrix:"<<"\n";
            printM(cmatrix); //Print new matrix
            store(cmatrix,getUnM(num)); //Store matrix
            system("PAUSE");
        }
        else if(val == 9) //Inverse Matrix
        {
            int num;
            cout<<"Which matrix do you want to invert?\n>";
            cin>>num;
            //Changes original/unmod matrix to matrix stored in num value.
            setMtx(num,getMtx(num));
            printUM(getUnM(num));
            cmatrix.resize(getMtx(num).numrows(), getMtx(num).numcols()); //resize to fit target matrix
            cmatrix = inverse(getMtx(num),getUnM(num)); //matrix = target matrix
            cout<<"Inverse Matrix:"<<"\n";
            printM(cmatrix); //Print new matrix
            store(cmatrix,getUnM(num)); //Store matrix
            system("PAUSE");
        }
        else if(val == 10) //Linear Regression
        {
            rVals();
            //printM(dataY);
            apmatrix<double>reg;
            reg.resize(1,2);
            reg = linReg(dataX, dataY);
            cout<<"y=ax+b => y="<<reg[0][1]<<"x+"<<reg[0][0];
            system("PAUSE");
        }
        else if(val == 11) //Matrix Regression
        {
            rVals();
            cout<<"Enter operations. Type 'm' for scalar multiplication, 'p' for power, 'sqrt' for square root, 'l' for log, 'sin',  'cos', 'tan', or 'd' when finished.\n";
            bool eF = true;
            int row=0;
            dataF.resize(row+1, 1);
            dataV.resize(row+1, 1);
            while(eF)
            {
                cout<<">";
                string n;
                cin>>n;
                if(n=="d")
                {
                    dataF.resize(dataV.numrows()-1,1);
                    dataV.resize(dataF.numrows(),1);
                    eF = false;
                }

                else
                {
                    if(n=="m")
                    {
                        cout<<"Enter scalar value: ";
                        cin>>dataV[row][0];
                    }
                    else if(n == "p")
                    {
                        cout<<"Enter exponent: ";
                        cin>>dataV[row][0];
                    }
                    else if(n == "l")
                    {
                        cout<<"Enter log base value: ";
                        cin>>dataV[row][0];
                    }
                    else
                    {
                        dataV[row][0] = 0;
                    }
                    double num;
                    // dataF[row][0]=strtod(n.c_str(),NULL);
                    dataF[row][0] = n;
                    row++;
                    dataF.resize(row+1,1);
                    dataV.resize(row+1,1);
                }
            }
            apmatrix<double> view = regrssn(dataX, dataY, dataF, dataV);
            /*for(int i = 0; i<view.numrows(); i++) {
                for(int j = 0; j<view.numcols(); j++) {
                    cout<<view[i][j];
                }
                cout<<endl;
            }*/
            cout<<"y="<<view[0][0]<<"+"<<view[1][0]<<"x+";
            for(int i = 2; i<view.numrows(); i++)
            {
                for(int j = 0; j<view.numcols(); j++)
                {
                    if(dataF[i-2][0]==("m") || dataF[i-2][0]==("M"))   //Scalar multiplication
                    {
                        cout<<view[i][0]<<"*x";
                    }
                    else if(dataF[i-2][0]==("p") || dataF[i-2][0]==("P"))   //Power
                    {
                        cout<<view[i][0]<<"*";
                        cout<<"x^"<<dataV[i-2][0];
                    }
                    else if(dataF[i-2][0]==("sqrt"))   //Square Root
                    {
                        cout<<view[i][0]<<"*sqrt(x)";
                    }
                    else if (dataF[i-2][0]==("l") || dataF[i-2][0]==("L"))   //Log Base
                    {
                        cout<<view[i][0]<<"*log("<<dataV[i-2][0]<<",x)";
                    }
                    else if (dataF[i-2][0]==("sin"))   //Sin
                    {
                        cout<<view[i][0]<<"*sin(x)";
                    }
                    else if (dataF[i-2][0]==("cos"))   //Cos
                    {
                        cout<<view[i][0]<<"*cos(x)";
                    }
                    else if (dataF[i-2][0]==("tan"))   //Tan
                    {
                        cout<<view[i][0]<<"*tan(x)";
                    }
                }
                if(i+2<view.numrows()&&view[i+1][0]>=0)
                {
                    cout<<"+";
                }
            }
        }
        else if(val == 12) //QR Factorization
        {
            int num;
            cout<<"Which matrix do you want to factorize?\n>";
            cin>>num;
            //Changes original/unmod matrix to matrix stored in num value.
            setMtx(num,getMtx(num));
            printUM(getUnM(num));
            apmatrix<double>Q = orthonormal(orthogonal(getUnM(num)));
            apmatrix<double>R = mult(transpose(Q,Q),getUnM(num));
            apmatrix<double>eVals = eV(getUnM(num),Q,R);
            eVals.resize(1,getUnM(num).numrows());
            cout<<"EigenValues: ";
            for(int z=0; z<eVals.numcols();z++) //Prints Eigenvalues
                cout<<eVals[0][z]<<" ";
                cout<<endl<<"Done. "<<spent<<" seconds.";
            system("PAUSE");
        }
        else if(val == 0) //Quit
        {
            // onClose();
            return EXIT_SUCCESS;
        }
        else
        {
            cout<<"Please enter a 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, or 0.\n";
        }
    }//I
    system("PAUSE");
    return EXIT_SUCCESS;
}
void enter()
{
    //Declare Matrix Size
    apmatrix<double>mtx;
    apmatrix<double>orig;
    cout<<"Enter the 'row' size of the matrix: "<<"\t";
    cin>>row;
    cout<<"Enter the 'col' size of the matrix: "<<"\t";
    cin>>col;
    //Enter Values for Matrix
    mtx.resize(row, col);
    orig.resize(row,col);
    cout<<"Press 'd' to go back one element"<<endl;
    for(int f=0; f<row; f++) //A
    {
        for(int j=0; j<col; j++) //1
        {
            string cord;
            cout<<"Enter value for ("<<(f+1)<<","<<(j+1)<<"): ";
            cin>>cord;
            if(cord == "d") {
                if(j==0){
                    f--;
                    j=col-2;
                }
                else
                    j-=2;
                continue;
            }
            double coord = strtod(cord.c_str(),NULL);
            mtx[f][j]=coord;
            orig[f][j]=coord;
        } //1
    } //A
    //Overwrite or Store
    bool st;
    while(st)
    {
        cout<<"Which matrix would you like to store this in (1-6)? Enter '0' to go to menu.\n>";
        int indx;
        cin>>indx;
        if(indx>=1 && indx<=6)
        {
            setMtx(indx,mtx);
            cout<<"\nThis matrix is labeled: "<<indx<<"\n";

            st = false;
        }
        else if(indx==0)
            st=false;
        else
        {
            cout<<"\nPlease enter a whole number value from 1-6\n";
        }
    }
}
void setMtx(int indx, apmatrix<double>mtx)
{
    if(indx == 1)
    {
        m1.resize(row,col);
        m1 = mtx;
        unmod1.resize(row,col);
        unmod1 = mtx;
    }
    else if(indx == 2)
    {
        m2.resize(row,col);
        m2=mtx;
        unmod2.resize(row,col);
        unmod2 = mtx;
    }
    else if(indx == 3)
    {
        m3.resize(row,col);
        m3 = mtx;
        unmod3.resize(row,col);
        unmod3 = mtx;
    }
    else if(indx == 4)
    {
        m4.resize(row,col);
        m4 = mtx;
        unmod4.resize(row,col);
        unmod4 = mtx;
    }
    else if(indx == 5)
    {
        m5.resize(row,col);
        m5 = mtx;
        unmod5.resize(row,col);
        unmod5 = mtx;
    }
    else if (indx== 6)
    {
        m6.resize(row,col);
        m6 = mtx;
        unmod6.resize(row,col);
        unmod6 = mtx;
    }
    else
    {
        cout<<"Something is not right.";
        system("EXIT_FAILURE");
    }
}
apmatrix<double> getMtx(int val)
{
    if (val ==1)
    {
        return m1;
    }
    else if(val ==2)
    {
        return m2;
    }
    else if (val ==3)
    {
        return m3;
    }
    else if(val==4)
    {
        return m4;
    }
    else if(val==5)
    {
        return m5;
    }
    else if(val=6)
    {
        return m6;
    }
}
apmatrix<double> getUnM(int val)
{
    if (val ==1)
    {
        return unmod1;
    }
    else if(val ==2)
    {
        return unmod2;
    }
    else if (val ==3)
    {
        return unmod3;
    }
    else if(val==4)
    {
        return unmod4;
    }
    else if(val==5)
    {
        return unmod5;
    }
    else if(val=6)
    {
        return unmod6;
    }
}
void store(apmatrix<double>matrix, apmatrix<double>orig)
{
    bool h =true;
    while(h)
    {
        cout<<"Which matrix would you like to store this in?\nEnter 'c' for current, '0' to go to menu, or enter 1-6.\n>";
        char stor;
        cin>>stor;
        if(stor == 'c' || stor=='C')   //Stores matrix in state that it's already in.
        {
            bool test = true;
            for(int i=1; i<=6; i++) //Iterates through each matrix
            {
                if(getMtx(i).numrows() == orig.numrows() && getMtx(i).numcols() == orig.numcols())   //Checks that dimensions are equal
                {
                    for(int k = 0; k<orig.numrows(); k++)  //Unable to test condition: if(getMatrix(i) == orig)
                    {
                        for(int j = 0; j<orig.numcols(); j++)
                        {
                            if(getMtx(i)[k][j]!=orig[k][j]) //Individually checks each index until one does not match
                                test=false;
                        }
                    }
                }
                if(test)   //If all indexes are same, overwrites matrix and exits loop
                {
                    setMtx(i,matrix);
                    h=false;
                }
            }
            if(!test)
            {
                cout<<"[WARNING] This matrix does not exist in any existing slot. Please choose another option.";
            }
        }
        else if(stor>=49 && stor<=54)     //Int value of the char '1' - '6' (Look at ASCII table)
        {
            setMtx(stor-49,matrix);
            h=false;
        }
        else if(stor=='0') //quit
        {
            h=false;
        }
        else
        {
            cout<<"Please enter 'c' for current, '0' to go to menu or enter 1-6\n";
        }
    }
}

apmatrix<double> transpose(apmatrix<double> mtrix, apmatrix<double>orig)    //orig = unmod matrix
{
    apmatrix<double>trans;
    trans.resize(mtrix.numcols(), mtrix.numrows());
    // matrix.resize(orig.numrows(),orig.numcols()); //If not a square, flips m with n
    for(int i = 0; i<mtrix.numrows(); i++)
    {
        for(int j=0; j<mtrix.numcols(); j++)
        {

            trans[j][i] = mtrix[i][j]; //Switching the variables
        }
    }
    return trans;
}
apmatrix<double> add(apmatrix<double>mtx, apmatrix<double> matrix)
{
    apmatrix<double> ph;
    if(mtx.numrows() == matrix.numrows() && mtx.numcols() == matrix.numcols())
    {
        for(int i = 0; i<matrix.numrows(); i++)
        {
            for(int j = 0; j<matrix.numcols(); j++)
            {
                ph.resize(matrix.numrows(), matrix.numcols());
                ph[i][j] = mtx[i][j] + matrix[i][j];
            }
        }
        return ph;

    }
    else
    {
        apmatrix<double>m;
        // m.resize(1,1);
        m[0][0]=1.0;
        return m;
    }
    // return matrix;
}
apmatrix<double> subtract(apmatrix<double>mtx, apmatrix<double> mtrix)
{
    apmatrix<double> ph;
    if(mtx.numrows() == mtrix.numrows() && mtx.numcols() == mtrix.numcols())
    {
        for(int i = 0; i<mtrix.numrows(); i++)
        {
            for(int j = 0; j<mtrix.numcols(); j++)
            {
                ph.resize(mtrix.numrows(), mtrix.numcols());
                ph[i][j] = mtx[i][j] - mtrix[i][j];
            }
        }
        return ph;

    }
    else
    {
        apmatrix<double>m;
        // m.resize(1,1);
        m[0][0]=1.0;
        return m;
    }
    // return matrix;
}
apmatrix<double> mult(apmatrix<double>mtx, apmatrix<double> mtrix)
{
    apmatrix<double>dh;
    if(mtx.numcols() == mtrix.numrows())   //Checks if col of mtx1 = row of mtx2
    {
        dh.resize(mtx.numrows(), mtrix.numcols()); //Sets size (mtx1.col, mtx2.row);
        // cout<<"rows:rows: "<<dh.numrows()<<" cols: "<<dh.numcols()<<"\n";
       for(int i = 0; i<dh.numrows(); i++)
        {
            for(int j = 0; j<dh.numcols(); j++)   //This loop is for value to be entered
            {
                //double index=0;
                dh[i][j]=0; //Fixes Memory issue when method reiterates
                for(int k=0; k<mtx.numcols(); k++)
                {
                    dh[i][j]+=mtx[i][k]*mtrix[k][j];
                }
                //dh[i][j] = index;
            }
        }
        return dh;

    }
    else
    {
        apmatrix<double>m;
        // m.resize(1,1);
        m[0][0]=1.0;
        return m;
    }
}
apmatrix<double> rE(apmatrix<double> matrix, apmatrix<double> unmd)
{
//!//Row Echelon Form
    /*Variables:*/
    //cout<<"Matrix after import: "; printM(matrix);
    for(int r=0; r<matrix.numrows()-1; r++) //B //Target Row
    {
        for(int i =0; i<matrix.numrows(); i++)   //Checks if all leading values are 0.
        {
            if(matrix[i][0]!=0)
            {
                i=matrix.numrows()-1;
                continue;
            }
            else if (i==matrix.numrows()-1 && matrix[i][0]==0) {
                r++; //Skips target row for reduce

               bool nonZ = true;
               double simp = 1; //While the element is 0, divides by 1
                for(int h = 0; h<matrix.numcols(); h++) { //Simplifies previous Target Row in case col of 0s.
                    if(matrix[r-1][h]!=0 && nonZ) {
                        simp = matrix[r-1][h];
                    nonZ = false;
                    }
                        matrix[r-1][h]*=(double)1/simp; //Simplifies all elements, first non-zero and beyond in specified row.
                }

            }
        }
        if(matrix[r][r]==0)  //1 //[Diagonal] Leading 0 Checker
        {
            bool t = true;
            while(t)  //a //Instead of break statement.
            {
                for(int x=r+1; x<matrix.numrows(); x++) //i //Switches Target row if leading coefficient is 0
                {
                    if(matrix[x][r]!=0)  //! //Main Diagonal
                    {
                        double temp[matrix.numcols()]; //=matrix[r];
                        for(int h=0; h<matrix.numcols(); h++) //z //Flips rows if leading value is 0
                        {
                            temp[h] = matrix[0][h];
                            matrix[r][h]=matrix[x][h];
                            matrix[x][h]=temp[h];
                        }//z
                        t=false;
                    }//!
                }//i
            }//a
        }//1
        // cout<<"Matrix after 0 checker: "; printM(matrix);
        for(int i=r+1; i<matrix.numrows(); i++) //b
        {
            for(int j=matrix.numcols()-1; j>=0; j--) //i  //Simplification
            {
                if(matrix[i][r]!=0)  //!
                {
                    //Creates leading 0
                    //   cout<<"Math: "<<matrix[i][j]<<"="<<matrix[r][j]<<"*"<<matrix[i][r]<<"/"<<matrix[r][r]<<"-"<<matrix[i][j];
                    matrix[i][j]=(double)((double)(matrix[r][j]*(double)(matrix[i][r]/matrix[r][r])*-1)+(double)matrix[i][j]); //Rt*(Rm/Rt)-Rm
                    // cout<<"\nr: "<<r<<"\ni: "<<i<<"\nj: "<<j;
                    // cout<<"="<<matrix[i][j]<<endl;
                }//!

            } //i
        } //b

    } //B
//printM(matrix);
    //Makes leading value 1
    for(int i = 0; i<matrix.numrows(); i++)
    {
        for(int j = matrix.numcols()-1; j>=i; j-- )
        {

            if(matrix[i][i]!=0)
            {
                //cout<<"Leading 1: "<<matrix[i][j]<<"="<<matrix[i][j]<<"/"<<matrix[i][i]<<"=";
                double temp = (double)(matrix[i][j]/matrix[i][i]); //Simplifies all elements in row by 1/first element
                matrix[i][j] = temp;
                //cout<<matrix[i][j]<<endl;
            }
        }
    }
    return matrix;
}
apmatrix<double> rrE(apmatrix<double> matrix, apmatrix<double> unmod)
{

    for(int r=0; r<matrix.numrows()-1; r++) //Targeted Row
    {
        for(int i=r+1; i<matrix.numrows(); i++) //Alter row
        {
            bool r0=false;
           for(int j = 0; j<matrix.numcols(); j++) { //Row of 0s checker
                if(matrix[i][j]!=0)
                    break;
                else if(j==matrix.numcols()-1)
                        r0=true;
           }
                if(r0)
                    continue;
            for(int j=r+1; j<matrix.numcols(); j++) //Col of Alter row
            {
                printM(matrix);
                //if(unmod[i][r+1]!=0) { //Makes sure leading value in targeted !=0
                matrix[r][j] -= unmod[i][j]*unmod[r][i]/unmod[i][i]; //B=B-E(Rt/Et)
                cout<<matrix[r][j]<<endl;
                //}
            }
            unmod = matrix; //Resets targeted row.
        }
    }
    return matrix;
}
apmatrix<double> inverse(apmatrix<double> mtx, apmatrix<double> unmd)
{
    apmatrix<double>invert;
    invert.resize(mtx.numrows(),mtx.numcols());
    mtx.resize(mtx.numrows(),(2*mtx.numcols())); //Space to add Identity
//Adds Identity to end of matrix
    for(int i = 0; i<mtx.numrows(); i++)
    {
        for(int j = (mtx.numcols()/2); j<mtx.numcols(); j++)
        {
            if(i==j-(mtx.numcols()/2))   //Finding where to put leading 1
            {
                mtx[i][j]=1;
                //  invert[i][j-mtx.numcols()/2]=1; //Becomes identity
            }
            else
            {
                mtx[i][j] = 0; //If not spot for leading 1, set index to 0.
                // invert[i][j-(mtx.numcols()/2)]=0;
            }
        }
    }
//matrix = mtx;
// mtx=rE(mtx,mtx);
    //cout<<"RE: "; printM(rE(mtx,mtx));
    mtx = rrE(rE(mtx,unmd),rE(mtx,unmd));//Reduces to get inverse on other side
    for(int i = 0; i<mtx.numrows(); i++)
    {
        for(int j = (mtx.numcols()/2); j<mtx.numcols(); j++)
        {
            invert[i][j-(mtx.numcols()/2)]=mtx[i][j]; //Sets invert to inverted matrix
        }
    }

    return invert;
}
//Prints matrix before operation
void printUM(apmatrix<double> unmod)
{
    int row = unmod.numrows();
    int col = unmod.numcols();

    cout<<"\n"<<"Original Matrix:"<<"\n"<<"[";
    //Formatting
    for(int i=0; i<row; i++) //D
    {
        for(int j=0; j<col; j++) //1
        {
            if(j==0&&i!=0)  //a
            {
                cout<<" ";
            }//a
            cout<<unmod[i][j];

            if(j<col-1) //b
            {
                cout<<" ";
            }//b
            else if (i<row-1)  //c
            {
                cout<<"\n";
            }//c
        }//2
    }//D
    cout<<"]";
    cout<<"\n"<<"\n";
}
//Prints matrix after operation
void printM(apmatrix<double> matrix)
{
    //Print out REF Matrix
    int row = matrix.numrows();
    int col = matrix.numcols();

    cout<<"[";
    for(int i=0; i<row; i++) //E
    {
        for(int j=0; j<col; j++) //1
        {
            if(j==0&&i!=0)  //a
            {
                cout<<" ";
            }//a
            cout<<matrix[i][j];

            if(j<col-1) //b
            {
                cout<<" ";
            }//b
            else if (i<row-1)  //c
            {
                cout<<"\n";
            }//c
        }//1
    }//E
    cout<<"]";
    cout<<"\n";
}
apmatrix<double> linReg(apmatrix<double>x, apmatrix<double>y)   //Linear Regression
{
    double a=0,b=0,sigX=0,sigY=0,sigX2=0, sigXY=0;
    int n = x.numrows(); //Sig = sigma
    for(int i = 0; i<x.numrows(); i++)
    {
        for(int j = 0; j<y.numcols(); j++)   //Only 1 column
        {
            sigX+=x[i][j]; //Sum of all x values
            sigX2+=pow(x[i][j],2); //Sum of all x^2 valued.
            sigY+=y[i][j]; //Sum of all y values
            sigXY+=x[i][j]*y[i][j]; //Sum of all x*y values
        }
    }
    b=(double)(sigY*sigX2-sigX*sigXY)/((double)n*sigX2-pow(sigX,2)); //formula for y-intercept
    a=(double)((double)n*sigXY-sigX*sigY)/((double)n*sigX2-pow(sigX,2)); //formula for slope
    apmatrix<double>test;
    test.resize(1,2); //Just for logistic purposes
    test[0][0]=a;
    test[0][1]=b;
    return test; //Inefficient
}
apmatrix<double> regrssn(apmatrix<double> x, apmatrix<double> y, apmatrix<string> f/*function*/, apmatrix<double> v /*degree of f*/)   //Z = [(MtM)^-1] * Mt * Y
{
    apmatrix<double>m; //Gets rid of overwrite bug.
    //x.resize(x.numrows(), x.numcols()+1);
    m.resize(x.numrows(), f.numrows()+2);

    for(int i = 0; i<x.numrows(); i++)   //Adds row of 1s to x matrix.
    {
        m[i][0] = 1;
        m[i][1] = x[i][0];
        for(int j = 2; j<m.numcols(); j++)   //Placeholder for uninitialized columns
        {
            m[i][j]=0;
        }
    }
    for(int i = 0; i<f.numrows(); i++)   //Number of functions to add to M matrix
    {
        //   x.resize(x.numrows(), x.numcols()+1); //Adds another column to x functions (M)
        for(int j = 0; j<x.numrows(); j++)
        {
            if(f[i][0]==("m") || f[i][0]==("M"))   //Scalar multiplication
            {
                m[j][i+2] =  v[i][0]*x[j][0];
            }
            else if(f[i][0]==("p") || f[i][0]==("P"))   //Power
            {
                m[j][i+2] = pow(x[j][0],v[i][0]);
            }
            else if(f[i][0]==("sqrt"))   //Square Root
            {
                m[j][i+2] = sqrt(x[j][0]);
            }
            else if (f[i][0]==("l") || f[i][0]==("L"))   //Log Base
            {
                m[j][i+2] = log(x[j][0])/log(v[i][0]);
            }
            else if (f[i][0]==("sin"))   //Sin
            {
                m[j][i+2] = sin(x[j][0]);
            }
            else if (f[i][0]==("cos"))   //Cos
            {
                m[j][i+2] = cos(x[j][0]);
            }
            else if (f[i][0]==("tan"))   //Tan
            {
                m[j][i+2] = tan(x[j][0]);
            }
            else
            {
                continue;
            }
        }
    }
    apmatrix<double> M = m;
    apmatrix<double> Mt = transpose(m,m);
    apmatrix<double> MtM = mult(Mt, M);
    y.resize(y.numrows()-1,y.numcols()); //Extra row of a 0 for some reason.
    apmatrix<double>MtY = mult(Mt,y);
    //apmatrix<double> MtMI = inverse(MTM,MTM);

    // cout<<"M: "<<M.numrows()<<"x"<<M.numcols()<<endl;
    // cout<<"Mt: "<<Mt.numrows()<<"x"<<Mt.numcols()<<endl;
    // cout<<"y: "<<y.numrows()<<"x"<<y.numcols()<<endl;
    // cout<<"MtM: "<<MtM.numrows()<<"x"<<MtM.numcols()<<endl;
    // cout<<"Mty: "<<MtY.numrows()<<"x"<<MtY.numcols()<<endl;
    return redaug(augment(MtM, MtY),MtM, MtY); //[MtM : MtY]

}
void rVals()   //Enter X and Y values for Matrix Regression
{
    cout<<"Enter x values. Type 'd' when finished.\n";
    bool eX = true;
    dataX.resize(row+1, 1);
    int row=0;
    while(eX)
    {
        cout<<">";

        string n;
        cin>>n;
        if(n=="d")
        {
            dataX.resize(dataX.numrows()-1,1);
            eX = false;
        }
        else
        {
            double num;
            dataX[row][0]=strtod(n.c_str(),NULL);
            row++;
            dataX.resize(row+1,1);
        }
    }
    //printM(dataX);
    cout<<"Enter y values.\n";
    dataY.resize(row+1, 1);
    row=0;
    for(int i = 0; i<dataX.numrows(); i++)
    {
        cout<<">";

        string n;
        cin>>n;
        if(n=="d")
        {
            dataY.resize(dataY.numrows()-1,1);

        }
        else
        {
            double num;
            dataY[row][0]=strtod(n.c_str(),NULL);
            row++;
            dataY.resize(row+1,1);
        }
    }
}
apmatrix<double> augment(apmatrix<double> A, apmatrix<double> B)   //Instead of inverse for Matrix Regression
{
    apmatrix<double>aug;

    aug = A;
//A = nxn
//B = nx1
    aug.resize(A.numrows(),A.numcols()+B.numcols()); //A augmented with B matrix
//A.resize(A.numrows(),(2*mtx.numcols())); //Space to add B to A
//Adds B to A
    for(int i = 0; i<A.numrows(); i++)
    {
        for(int j = A.numcols(); j<aug.numcols(); j++)
        {
            aug[i][j]=B[i][j-A.numcols()]; //Adds B to aug
        }

    }
    return aug;
}
apmatrix<double> redaug(apmatrix<double> aug, apmatrix<double> A, apmatrix<double> B)   //Reduces Augmented Matrix and gets new inverted matrix on right side
{
    apmatrix<double>X;
    X.resize(B.numrows(), B.numcols());
    aug = rrE(rE(aug,aug),rE(aug,aug));//Reduces to get X  matrix on right side
    for(int i = 0; i<aug.numrows(); i++)
    {
        for(int j = (A.numcols()); j<aug.numcols(); j++)
        {
            X[i][j-(A.numcols())]=aug[i][j]; //Sets X to x values
        }
    }
    return X;
}
apmatrix<double> proj(apmatrix<double> v, apmatrix<double> u)   //Projection of u onto v
{
//proj(v,u) = [(u dot v)/(v dot v)]*v
    double a=0,b=0,c=0; //Fixes memory issue
    a = mult(transpose(u,u),v)[0][0]; b = mult(transpose(v,v),v)[0][0]; c = (a/b);
    apmatrix<double>k;
    k.resize(v.numrows(),v.numcols());
    for(int i = 0; i<v.numrows(); i++)
    {
        for(int j = 0; j<v.numcols(); j++)
        {
            k[i][j]=v[i][j]*c;
        }
    }
    return k;
}
apmatrix<double> orthogonal(apmatrix<double> A)   //Orthogonal
{
    /*
    n = (number of columns of A)-1
    x = columns of A
    v0 = x0
    v1 = x1 - proj(v0,x1)
    vn = xn - proj(v0,xn) - proj(v1,xn) - ... - proj(v(n-1),xn)
    V = {v1, v2, ..., vn} or [v0 v1 ... vn]
    */
    apmatrix<double> V, x;
    int n = A.numcols();
    V.resize(A.numrows(),n);
    x.resize(A.numrows(), 1);
    for(int i = 0; i<A.numrows(); i++)
    {
        x[i][0]=A[i][1];
        V[i][0]=A[i][0];
    }
    for (int c = 1; c<n; c++)   //Iterates through each col of A as if each was its own matrix
    {
        apmatrix<double>vn,vc; //vn = Orthogonalized v (avoiding matrix overwriting of v); vc = previously orthogonalized v
       // vn.resize(v.numrows(), 1);
       // for(int l = 0; l<vn.numrows(); l++)
            vn=x;
        vc.resize(x.numrows(), 1);
        for(int i=0; i<c; i++)   //Vn = an-(sigma(t=1, n-1, proj(vt, xn))
        {
            for(int k = 0; k<V.numrows(); k++)
                vc[k][0] = V[k][i]; //Sets vc to designated v matrix
          //  cout<<"c: "<<c<<"\ni: "<<i<<"\nvc x proj\n";
            apmatrix<double>temp = proj(vc, x);
            for(int j = 0; j<A.numrows(); j++)
            {
               // cout<<vc[j][0]<<" "<<x[j][0]<<" "<<temp[j][0]<<endl;
                //vn = add(v,proj(vc,x));
                vn[j][0]-=temp[j][0]; //orthogonalize matrix
            }
           // cout<<endl;
        }
        for(int k = 0; k<V.numrows(); k++)
        {
            V[k][c]=vn[k][0]; //Subtracts orthogonalized col to V
        }
        if((c+1)<n)  //Matrix Out of Bounds Checker
        {
            for(int k = 0; k<A.numrows(); k++)
            {
                vn[k][0]=0;
                vc[k][0]=0;
                x[k][0]=A[k][c+1]; //Moves x onto next v
            }
        }
    }
    return V;
}
apmatrix<double> orthonormal(apmatrix<double> V)   //Orthonormal
{
// v/||v||
    apmatrix<double> O = V;
    for(int c = 0; c<V.numcols(); c++)  //Iterates through each column of V
    {
        double mag=0; //Magnitude of v
        for(int i = 0; i<O.numrows(); i++)
        {
            mag+=pow(O[i][c],2); //sqrt(v1n^2 + v2n^2 + ... + vmn^2)
        }
        for(int i=0; i<O.numrows(); i++)   //V/||v||
        {
            O[i][c]=V[i][c]/sqrt(mag);
        }
    }
    return O;
}
apmatrix<double>eV(apmatrix<double>A, apmatrix<double> Q, apmatrix<double>R1) {
int n = A.numrows(), m=A.numcols();
apmatrix<double>A1=A,A2,A3;
chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
for(int i = 1; i<(500*m); i++) { //Am=QtA(m-1)Q
       if(i!=1) {
        for(int s = 0; s<n; s++) {
            for(int d = 0; d<m; d++) {
                A3[s][d]=0;
            }
        }
        }
        if(i%2==1) {
        A2=mult(mult(transpose(Q,Q),A1),Q); //A2 = QtA1Q
        Q=orthonormal(orthogonal(A2)); //Q = orthonormal A2;
        A3=A2; //Used for final value

         for(int s = 0; s<A1.numrows(); s++) {
            for(int d = 0; d<A1.numcols(); d++) {
                A1[s][d]=0;
            }
        }

        }
        else {
            A1=mult(mult(transpose(Q,Q),A2),Q); //A1 = QtA2Q
            Q=orthonormal(orthogonal(A1)); //Q = orthonormal A1
            A3=A1; //Used for final value

             for(int s = 0; s<A2.numrows(); s++) {
            for(int d = 0; d<A2.numcols(); d++) {
                A2[s][d]=0;
            }
        }

        }
       /*cout<<"Matrix R:"<<endl;
        printM(mult(transpose(Q,Q),A2));
        cout<<"Matrix Q: "<<endl;
        printM(Q);
        cout<<endl;
        cout<<i<<":"<<endl;
        printM(A3);*/
}

chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
spent = (double)(chrono::duration_cast<chrono::microseconds>(end - begin).count())/1000000.0;
//spent = double(chrono::duration_cast<chrono::nanoseconds> (end - begin).count())/1000000.0;
//printM(A3);
apmatrix<double> eVals;
eVals.resize(1,m);
for(int i = 0; i<m; i++) { //Getting EVs on diagonal
eVals[0][i]=A3[i][i];
}
return eVals;
}
/*void onOpen() { //Imports Matrix from text file
ifstream file;
file.open("matrices.txt");
stringstream buffer;
buffer<<file.rdbuf(); //Allocates mem in file to buffer
string txt = buffer.str(); //Turns contents of buffer into single line string


size_t pos1=0, pos2; //Parameters to read

    bool r1 = true;
    while(r1) { //Iterates through first row

        int i=0;
        bool c = true; //Iterates through col
        while(c) {
              system("PAUSE");
            int j = 0;
            pos2 = txt.find(",",pos1); //Splits col indexes with ','
           // cout<<txt.substr(pos1, (pos2-pos1));
           //double a = atof(txt.substr(pos1, (pos2-pos1)));
           string a = txt.substr(pos1, (pos2-pos1));
           cout<<a;
          // string loop;
         //  cin>>loop;
           //dataX[i][j] = strtod(b,NULL);
           pos1 = pos2+1;
            if(txt.substr((pos1+2),(pos1+2))==">") { //Designates new matrix
                pos1+=3; //Skips new line arg goes to next line
                r1=false;
                c=false;
            }
            else if(txt.substr((pos1+1),(pos1+1))==";") {//Designates new row
                pos1+=2;
                c=false;

            }
        j++;}
    i++;}

bool r2;
    while(r2) { //iterates through second row
        int i=0;
        bool c = true;
        while(c) {
            int j = 0;
            pos2 = txt.find(",",pos1);
            //m2[i][j] = atof(txt.substr(pos1, (pos2-pos1)));
            string a = txt.substr(pos1, (pos2-pos1));
            //cout<<"\n"<<a;
          //  dataY = strtod(a,NULL);
           pos1 = pos2+1;
            if(txt.substr((pos1+2),(pos1+2))==">") {
                pos1+=3;
                r1=false;
                c=false;
            }
            else if(txt.substr((pos1+1),(pos1+1))==";") {
                pos1+=2;
                c=false;
            }
        j++;}
    i++;}

    /*bool r3;
    while(r3) {
        int i=0;
        bool c = true;
        while(c) {
            int j = 0;
            pos2 = txt.find(",",pos1);
           // m3[i][j] =  atof(txt.substr(pos1, (pos2-pos1)));
           pos1 = pos2+1;
            if(txt.substr((pos1+2),(pos1+2))==">") {
                pos1+=3;
                r1=false;
                c=false;
            }
            else if(txt.substr((pos1+1),(pos1+1))==";") {
                pos1+=2;
                c=false;
            }
        j++;}
    i++;}
    bool r4;
    while(r4) {
        int i=0;
        bool c = true;
        while(c) {
            int j = 0;
            pos2 = txt.find(",",pos1);
           // m4[i][j] = atof(txt.substr(pos1, (pos2-pos1)));
           pos1 = pos2+1;
            if(txt.substr((pos1+2),(pos1+2))==">") {
                pos1+=3;
                r1=false;
                c=false;
            }
            else if(txt.substr((pos1+1),(pos1+1))==";") {
                pos1+=2;
                c=false;
            }
        j++;}
    i++;}
    bool r5;
    while(r5) {
        int i=0;
        bool c = true;
        while(c) {
            int j = 0;
            pos2 = txt.find(",",pos1);
            //m5[i][j] =  atof(txt.substr(pos1, (pos2-pos1)));
           pos1 = pos2+1;
            if(txt.substr((pos1+2),(pos1+2))==">") {
                pos1+=3;
                r1=false;
                c=false;
            }
            else if(txt.substr((pos1+1),(pos1+1))==";") {
                pos1+=2;
                c=false;
            }
        j++;}
    i++;}
    bool r6;
    while(r6) {
        int i=0;
        bool c = true;
        while(c) {
            int j = 0;
            pos2 = txt.find(",",pos1);
            //m6[i][j] =  atof(txt.substr(pos1, (pos2-pos1)));
           pos1 = pos2+1;
            if(txt.substr((pos1+2),(pos1+2))==">") {
                pos1+=3;
                r1=false;
                c=false;
            }
            else if(txt.substr((pos1+1),(pos1+1))==";") {
                pos1+=2;
                c=false;
            }
        j++;}
    i++;}
    file.close();
}
*/
/*
void onClose() {
    ofstream file;
    file.open("matrices.txt");
   // file.seekp(4);
     for(int i=0; i<row;i++) {
        for(int j=0;j<col;j++) {
            file<<m1[i][j]<<",";
        }
        file<<";";
       }
       file<<">";
       cout<<endl;

        for(int i=0; i<row;i++) {
        for(int j=0;j<col;j++) {
            file<<m2[i][j]<<",";
        }
        file<<";";
       }
        file<<">";
       cout<<endl;

        for(int i=0; i<row;i++) {
        for(int j=0;j<col;j++) {
            file<<m3[i][j]<<",";
        }
        file<<";";
       }
        file<<">";
       cout<<endl;
        for(int i=0; i<row;i++) {
        for(int j=0;j<col;j++) {
            file<<m4[i][j]<<",";
        }
        file<<";";
       }
        file<<">";
       cout<<endl;

        for(int i=0; i<row;i++) {
        for(int j=0;j<col;j++) {
            file<<m5[i][j]<<",";
        }
        file<<";";
       }
        file<<">";
       cout<<endl;
        for(int i=0; i<row;i++) {
        for(int j=0;j<col;j++) {
            file<<m6[i][j]<<",";
        }
        file<<";";
       }
        file<<">";

   file.close();
}
*/
