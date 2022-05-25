//#include<iostream>
//#include<cmath>
/*
Last modified on jan 4, 2016;
Changes:
MAJOR BUG FIX. Added initialization to 0 in the constructor.

oct 2, 2015.
Changes:
Overloaded Matrix * / scalar operations
Included gaussEliminationMatrix
*/

class complex
{
	public:
	double re,im;

	complex(double a=0.0,double b=0.0)
	{re=a;im=b;};

	double abs(void)
	{return sqrt(pow(re,2)+pow(im,2));}

	double ang(void)
	{return atan2(im,re);}

	complex conj(void)
	{
		complex result;
		result.re=re;
		result.im=-im;
		return result;
	}
};

complex operator+(complex first,complex second)
{
	complex result;
	result.re=first.re+second.re;
	result.im=first.im+second.im;
	return result;
};

complex operator-(complex first,complex second)
{
	complex result;
	result.re=first.re-second.re;
	result.im=first.im-second.im;
	return result;
};

complex operator*(complex first,complex second)
{
	complex result;
	result.re=first.re*second.re-first.im*second.im;
	result.im=first.re*second.im+first.im*second.re;
	return result;
};

complex operator/(complex first,complex second)
{
	complex num;
	num=first*second.conj();
	num.re=num.re/pow(second.abs(),2);
	num.im=num.im/pow(second.abs(),2);
	return num;
};


std::ostream& operator<<(std::ostream& output,const complex &num)
{
	return output<<num.re<<","<<num.im<<"i ";
}


template <typename T, int N> 
class arr
{
	public:
	T myarray[N];

	int max;
	arr()
	{
		max=N;
		for (int i=0;i<N;i++) myarray[i]= (T) 0;
	};

	T operator[](const int i) const {return myarray[i];}

	T & operator[](const int i)
	{return myarray[i];}

	arr<T,N> &operator= (const T &rhs)
	{
		for (int i=0;i<N;i++)
		{myarray[i]=rhs;};
		return *this;
	};

};

template<typename T,int N>
std::ostream& operator<<(std::ostream& output,const arr<T,N> &arrp)
{
	for (int i=0;i<N;i++)
	{
		output<<arrp[i]<<" ";
	};
	return output;
};

template <typename T,int N>
arr<T,N> operator+(arr<T,N> first,arr<T,N> second)
{
arr<T,N> result;
for (int i=0;i<N;i++)
{result[i]=first[i]+second[i];};
return result;
};

template <typename T,int N>
arr<T,N> operator-(arr<T,N> first,arr<T,N> second)
{
	arr<T,N> result;
	for (int i=0;i<N;i++)
	{result[i]=first[i]-second[i];};
	return result;
};

template <typename T,int N>
arr<T,N> operator*(arr<T,N> first,arr<T,N> second)
{
	arr<T,N> result;
	for (int i=0;i<N;i++)
	{result[i]=first[i]*second[i];};
	return result;
};

template <typename T,int N>
arr<T,N> operator*(arr<T,N> first,T scalar)
{
	arr<T,N> result;
	for (int i=0;i<N;i++)
	{result[i]=first[i]*scalar;};
	return result;
};

template <typename T,int N>
arr<T,N> operator*(T scalar,arr<T,N> first)
{
	arr<T,N> result;
	for (int i=0;i<N;i++)
	{result[i]=first[i]*scalar;};
	return result;
};

template <typename T,int N>
arr<T,N> operator/(arr<T,N> first,T scalar)
{
	arr<T,N> result;
	for (int i=0;i<N;i++)
	{result[i]=first[i]/scalar;};
	return result;
};

template <typename T,int N>
arr<T,N> operator/(arr<T,N> first,arr<T,N> second)
{
	arr<T,N> result;
	for (int i=0;i<N;i++)
	{result[i]=first[i]/second[i];};
	return result;
};

template <typename T, int N1, int N2>
class mat
{
	public:
	arr <T,N2> mymat[N1];
	int imax,jmax;

	mat()
	{
		imax=N1;
		jmax=N2;

		for (int i=0;i<N1;i++)
		{
			for (int j=0;j<N2;j++)
			{
				mymat[i][j]= (T) 0;
			}
		}
	};

	bool symmetric(void)
	{
		bool symm;
		for (int i=0;i<N1;i++)
		{
			for (int j=0;j<i;j++)
			{
				if (mymat[i][j]==mymat[j][i])
				{
					symm=true;
				}
				else
				{
					symm=false;
					goto exitline;
				};
			};
		};
		exitline:return symm;
	};


	T operator()(const int i,const int j) const {return mymat[i][j];}

	T & operator()(const int i,const int j)
	{return mymat[i][j];}

	mat<T,N1,N2> &operator= (const T &rhs)
	{
		for (int i=0;i<N1;i++)
		{
			for (int j=0;j<N2;j++)
			{mymat[i][j]=rhs;};
		};
	return *this;
	};

};

template<typename T,int N1, int N2>
std::ostream& operator<<(std::ostream& output,const mat<T,N1,N2> &matp)
{
	for (int i=0;i<N1;i++)
	{
		for (int j=0;j<N2;j++)
		{output<<matp(i,j)<<" ";};
		output<<std::endl;
	};
	return output;
};


template <typename T,int N1,int N2>
mat<T,N1,N2> operator+(mat<T,N1,N2> first,mat<T,N1,N2> second)
{
	mat<T,N1,N2> result;
	for (int i=0;i<N1;i++)
	{
		for (int j=0;j<N2;j++)
		{
			result(i,j)=first(i,j)+second(i,j);
		};
	};
	return result;
};

template <typename T,int N1,int N2>
mat<T,N1,N2> operator-(mat<T,N1,N2> first,mat<T,N1,N2> second)
{
	mat<T,N1,N2> result;
	for (int i=0;i<N1;i++)
	{
		for (int j=0;j<N2;j++)
		{
			result(i,j)=first(i,j)-second(i,j);
		};
	};
	return result;
};

template <typename T,int N1,int N2>
mat<T,N1,N2> operator*(mat<T,N1,N2> first,mat<T,N1,N2> second)
{
	mat<T,N1,N2> result;
	for (int i=0;i<N1;i++)
	{
		for (int j=0;j<N2;j++)
		{
			result(i,j)=first(i,j)*second(i,j);
		};
	};
	return result;
};

template <typename T,int N1,int N2>
mat<T,N1,N2> operator/(mat<T,N1,N2> first,mat<T,N1,N2> second)
{
	mat<T,N1,N2> result;
	for (int i=0;i<N1;i++)
	{
		for (int j=0;j<N2;j++)
		{
			result(i,j)=first(i,j)/second(i,j);
		};
	};
return result;
};

template <typename T,int N1,int N2>
mat<T,N1,N2> operator*(mat<T,N1,N2> first,T scalar)
{
	mat<T,N1,N2> result;
	for (int i=0;i<N1;i++)
	{
		for (int j=0;j<N2;j++)
		{
			result(i,j)=first(i,j)*scalar;
		};
	};
	return result;
};

template <typename T,int N1,int N2>
mat<T,N1,N2> operator*(T scalar,mat<T,N1,N2> first)
{
	mat<T,N1,N2> result;
	for (int i=0;i<N1;i++)
	{
		for (int j=0;j<N2;j++)
		{
			result(i,j)=first(i,j)*scalar;
		};
	};
return result;
};

template <typename T,int N1,int N2>
mat<T,N1,N2> operator/(mat<T,N1,N2> first,T scalar)
{
	mat<T,N1,N2> result;
	for (int i=0;i<N1;i++)
	{
		for (int j=0;j<N2;j++)
		{
			result(i,j)=first(i,j)/scalar;
		};
	};
	return result;
};

template <typename T,int N1,int N2,int N3,int N4>
mat<T,N1,N4> matMul(mat<T,N1,N2> first,mat<T,N3,N4> second)
{

	mat<T,N1,N4> result;
	for (int i=0;i<N1;i++)
	{
		for (int j=0;j<N4;j++)
		{
			result(i,j)=0.0;
			for (int k=0;k<N2;k++)
			{
				result(i,j)=result(i,j)+first(i,k)*second(k,j);
			};
		};
	};
	return result;
};



template <int N>
class Factorial
{
	public:
	enum {value=N*Factorial<N-1>::value};
};

template <>
class Factorial<0>
{
	public:
	enum {value=1};
};

template<class T,int N>
inline T dot(arr<T,N> &a, arr<T,N> &b);

template<int m>
class metaDot
{
	public:
	template<class T,int N>
	static T f(arr<T,N> & a,arr<T,N> & b)
	{
		return a[m]*b[m]+metaDot<m-1>::f(a,b);
	}
};

template<>
class metaDot<0>
{
	public:
	template<class T,int N>
	static T f(arr<T,N> & a,arr<T,N> & b)
	{
		return a[0]*b[0];
	}
};

template<class T,int N>
inline T dot(arr<T,N> &a, arr<T,N> &b)
{
	return metaDot<N-1>::f(a,b);
};

template<class T,int N1,int N2>
mat <T,N2,N1> transpose(mat <T,N1,N2> A)
{
	mat <T,N2,N1> At;
	for (int i=0;i<N1;i++)
	{
		for (int j=0;j<N2;j++)
		{
			At(j,i)=A(i,j);
		};
	};
	return At;
}

template<class T,int N1,int N2,int N3, int N4>
mat<T,N3,N4> powerMethod(mat <T,N1,N2> A,mat <T,N3,N4> x, int N)
{
	T norm;
	for (int i=0;i<N;i++)
	{
		x=matMul(A,x);
		norm=sqrt(matMul(transpose(x),x)(0,0));
		for (int j=0;j<N3;j++)
		{
			x(j,0)/=norm;
		};
	};
	return x;
};

template<class T,int N1,int N2,int N3,int N4>
mat<T,N3,N4> gaussEliminationMatrix(mat <T,N1,N2> A, mat <T,N3,N4> B)
{
	mat <T,N3,1> b;
	mat <T,N3,1> x;
	mat <T,N3,N4> Result;

	for (int j=0;j<N4;j++)
	{

		for (int i=0;i<N3;i++)
		{
		  b(i,0)=B(i,j);
		};
		 
		x=gaussElimination(A,b);
		 
		for (int i=0;i<N3;i++)
		{
		  Result(i,j)=x(i,0);
		};
	 
	};
	return Result;
};

template<class T,int N1,int N2,int N3,int N4>
mat<T,N3,N4> gaussElimination(mat <T,N1,N2> A,mat <T,N3,N4> b)
{
	mat <T,N3,N4> x;
	mat <T,1,N2> matTemp;
	T rTemp,lTemp,multiplier= (T) 0;

	for (int row=0;row<N1;row++)
	{
		for (int k=0;k<N1;k++)
		{
			if (A(k,k)== (T) 0)
			{
				for (int i=0;i<N1;i++)
				{
					if (A(i,k)!= (T) 0)
					{
						{
							for (int l=0;l<N2;l++)
							{
								matTemp(0,l)=A(i,l);
								A(i,l)=A(k,l);
								A(k,l)=matTemp(0,l);
								rTemp=b(i,0);
								b(i,0)=b(k,0);
								b(k,0)=rTemp;
							};
							break;
						};
					};

				}; 
			};
		};
	};


	for (int j=0;j<N2;j++)
	{
		for (int i=j+1;i<N1;i++)
		{
			if (A(i,j)!=(T) 0)
			{
				multiplier=A(i,j)/A(j,j);
				for(int k=j;k<N2;k++)
				{
					A(i,k)=A(i,k)-A(j,k)*multiplier;
				};
				b(i,0)=b(i,0)-b(j,0)*multiplier;
			};
		};
	};


	for(int i=N1-1;i>=0;i--)
	{
		for (int j=i+1;j<N2;j++)
		{
			lTemp=lTemp+A(i,j)*x(j,0);
		};
		x(i,0)=(b(i,0)-lTemp)/A(i,i);
		lTemp=(T) 0;
	};

	return x;
};

template<class T,int N1,int N2>
void LUfactorization(mat <T,N1,N2> A,mat <T,N1,N2> &L,mat <T,N1,N2> &U)
{
	T multiplier= (T) 0;
	L= (T) 0;
	U= (T) 0;

	for (int i=0;i<N1;i++)
	{
		L(i,i)= (T) 1;
	};
	for (int j=0;j<N2;j++)
	{
		for (int i=j+1;i<N1;i++)
		{
			if (A(i,j)!=(T) 0)
			{
				multiplier=A(i,j)/A(j,j);
				for(int k=j;k<N2;k++)
				{
					A(i,k)=A(i,k)-A(j,k)*multiplier;
				};
				L(i,j)=multiplier;
			};
		};
	};
	U=A;
};

template <class T,int N1,int N2,int N3,int N4>
mat <T,N3,N4> Cholesky(mat <T,N1,N2> A, mat <T,N3,N4> b)
{
	mat <T,N1,N2> L,U;
	mat <T,N3,N4> x,y;
	T lTemp= (T) 0;

	LUfactorization(A,L,U);

	for(int i=0;i<N1;i++)
	{
		for (int j=0;j<i;j++)
		{
			lTemp=lTemp+L(i,j)*y(j,0);
		};
		y(i,0)=(b(i,0)-lTemp)/L(i,i);
		lTemp=(T) 0;
	};

	for(int i=N1-1;i>=0;i--)
	{
		for (int j=i+1;j<N2;j++)
		{
			lTemp=lTemp+U(i,j)*x(j,0);
		};
		x(i,0)=(y(i,0)-lTemp)/U(i,i);
		lTemp=(T) 0;
	};
	return x;
};

template <class T,int N1,int N2,int N3,int N4>
mat <T,N3,N4> Seidel(mat <T,N1,N2> A, mat <T,N3,N4> b,int N)
{
	mat <T,N3,N4> x;
	x= (T) 0;
	T lTemp= (T) 0;

	for (int iter=0;iter<N;iter++)
	{
		for (int i=0; i<N1; i++)
		{
			for (int j=0;j<N2;j++)
			{
				if (j!=i)
				{
					lTemp=lTemp+x(j,0)*A(i,j);
				}
			};
			x(i,0)=(b(i,0)-lTemp)/A(i,i);
			lTemp= (T) 0;
		};
	};
	return x;
};

template <class T,int N1,int N2,int N3,int N4>
mat <T,N3,N4> Jacobi(mat <T,N1,N2> A, mat <T,N3,N4> b,int N)
{
	mat <T,N3,N4> x,xTemp;
	x= (T) 0;
	T lTemp= (T) 0;

	for (int iter=0;iter<N;iter++)
	{
		for (int i=0; i<N1; i++)
		{
			for (int j=0;j<N2;j++)
			{
				if (j!=i)
				{
					lTemp=lTemp+x(j,0)*A(i,j);
				}
			};
			xTemp(i,0)=(b(i,0)-lTemp)/A(i,i);
			lTemp= (T) 0;
		};
		x=xTemp;
	};
	return x;
};
 


