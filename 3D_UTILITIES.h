#pragma once

#include <iostream>


#define rad_to_DEG(x) (180/M_PI*x)





	template <typename T,int size_rows,int size_columns>
	struct s_matrix
	{

		T data[size_rows * size_columns]={0};

		
		template <int size_rows2,int size_columns2>
		s_matrix<T,size_rows,size_columns2> operator*(s_matrix<T,size_rows2, size_columns2> mat2)
		{
			if (size_columns != size_rows2)
			{
				std::cerr << "The number of columns in the first matrix should be equal to the number of rows in the second : " << size_columns << " != " << size_rows2 << "\n";
				return s_matrix<T, size_rows, size_columns2>();
			}
			s_matrix<T, size_rows, size_columns2> newmat = { 0 };
			for (int i = 0; i < size_rows; i++)
				for (int k = 0; k < size_rows2; k++)
					for (int j = 0; j < size_columns2; j++)
						newmat.data[i * size_columns2 + j] += this->data[i * size_columns + k] * mat2.data[k * size_columns2 + j];
			return newmat;
		 }

		friend
			std::ostream& operator<<(std::ostream& os, s_matrix<T, size_rows,size_columns> mat)
		{
			
			for (int i = 0; i < size_rows; i++)
			{
				for (int j = 0; j < size_columns; j++)
					os << mat.data[i * size_columns + j] << " ";
				os << "\n";
			}
			
			return os;
		}

		s_matrix<T, size_rows, size_columns> operator/(T div)
		{
			s_matrix<T, size_rows, size_columns> cal;
			for (int i = 0; i < size_rows; i++)
				for (int j = 0; j < size_columns; j++)
					cal.data[i * size_columns + j] = data[i * size_columns + j] / div;
			return cal;
		}

		s_matrix<T, size_rows, size_columns> operator-(s_matrix<T, size_rows, size_columns> mat2)
		{
			s_matrix<T, size_rows, size_columns> cal;
			for (int i = 0; i < size_rows; i++)
				for (int j = 0; j < size_columns; j++)
					cal.data[i * size_columns + j]=	data[i * size_columns + j] - mat2.data[i * size_columns + j];
			return cal;
		}

		
		s_matrix<T, size_rows, size_columns> operator+(s_matrix<T, size_rows, size_columns> mat2)
		{
			s_matrix<T, size_rows, size_columns> cal;
			for (int i = 0; i < size_rows; i++)
				for (int j = 0; j < size_columns; j++)
					cal.data[i * size_columns + j] = data[i * size_columns + j] + mat2.data[i * size_columns + j];
			return cal;
		}

		s_matrix<T, size_rows, size_columns> operator-()
		{
			s_matrix<T, size_rows, size_columns> cal;
			for (int i = 0; i < size_rows; i++)
				for (int j = 0; j < size_columns; j++)
					cal.data[i * size_columns + j] = -data[i * size_columns + j] ;
			return cal;
		}

		void operator+=(s_matrix<T, size_rows, size_columns> mat2)
		{
			
			for (int i = 0; i < size_rows; i++)
				for (int j = 0; j < size_columns; j++)
					data[i * size_columns + j] += mat2.data[i * size_columns + j];
			
		}

		void operator-=(s_matrix<T, size_rows, size_columns> mat2)
		{

			for (int i = 0; i < size_rows; i++)
				for (int j = 0; j < size_columns; j++)
					data[i * size_columns + j] -= mat2.data[i * size_columns + j];
			
		}
	
		s_matrix<T, size_rows, size_columns> operator*(const float f)
		{
			s_matrix<T, size_rows, size_columns> cal;
			for (int i = 0; i < size_rows; i++)
				for (int j = 0; j < size_columns; j++)
					cal.data[i * size_columns + j] = data[i * size_columns + j]*f;
			return cal;

		}

		void operator*=(const float f)
		{
			
			for (int i = 0; i < size_rows; i++)
				for (int j = 0; j < size_columns; j++)
					data[i * size_columns + j] = data[i * size_columns + j] * f;
			
		}

	};


	typedef s_matrix<int, 4,4> imat4x4;
	typedef s_matrix<float, 4,4> fmat4x4;
	typedef s_matrix<int, 3,3> imat3x3;
	typedef s_matrix<float, 3,3> fmat3x3;

	struct m_4xRot : public fmat4x4
	{


		m_4xRot(float angle)
		{
			float sinx = sinf(angle);
			float cosx = cosf(angle);
			data[0] = 1;
			data[15] = 1;
			data[5] = cosx;
			data[6] = -sinx;
			data[9] = sinx;
			data[10] = cosx;

		}

	};


	struct m_4yRot : public fmat4x4
	{


		m_4yRot(float angle)
		{
			float sinx = sinf(angle);
			float cosx = cosf(angle);
			data[5] = 1;
			data[15] = 1;
			data[0] = cosx;
			data[2] = sinx;
			data[8] = -sinx;
			data[10] = cosx;

		}

	};


	struct m_4zRot : public fmat4x4
	{


		m_4zRot(float angle)
		{
			float sinx = sinf(angle);
			float cosx = cosf(angle);
			data[10] = 1;
			data[15] = 1;
			data[0] = cosx;
			data[1] = -sinx;
			data[4] = sinx;
			data[5] = cosx;

		}

	};

	struct m_3xRot: public fmat3x3 {

		m_3xRot(float angle)
		{
			float sinx = sinf(angle);
			float cosx = cosf(angle);
			data[0] = 1;
			data[4] = cosx;
			data[5] = -sinx;
			data[7] = sinx;
			data[8] = cosx;
		}


	};

	struct m_3yRot : public fmat3x3 {

		m_3yRot(float angle)
		{
			float sinx = sinf(angle);
			float cosx = cosf(angle);
			data[4] = 1;
			data[0] = cosx;
			data[2] = sinx;
			data[6] = -sinx;
			data[8] = cosx;
		}


	};

	struct m_3zRot : public fmat3x3 {

		m_3zRot(float angle)
		{
			float sinx = sinf(angle);
			float cosx = cosf(angle);
			data[8] = 1;
			data[0] = cosx;
			data[1] = -sinx;
			data[3] = sinx;
			data[4] = cosx;
		}


	};



	template<typename T>
	struct vector3d : public s_matrix<T,3,1>
	{

		vector3d()
		{
			this->data[0] = NULL;
			this->data[1] = NULL;
			this->data[2] = NULL;
			
		}

		vector3d(s_matrix<float, 3, 1> mat)
		{
			 
			for (int i = 0; i < 3; i++)
				this->data[i] = mat.data[i];
			

		}


		vector3d(T _X, T _Y, T _Z)
		{
			this->data[0] = _X;
			this->data[1] = _Y;
			this->data[2] = _Z;
			
		}
		
		float Dot_product(vector3d<T> Vector)
		{	
			
			return this->data[0]*Vector.data[0]+ this->data[1] * Vector.data[1]+ this->data[2] * Vector.data[2];
		}

		float Magnitude()
		{
		
			return sqrt(Dot_product(*this));
		}


		float angle_To(vector3d<T> Vector)
		{
			return acosf(this->Dot_product(Vector) / (this->Magnitude() * Vector.Magnitude()));

		}

		vector3d<T> cross_product(vector3d<T> Vector)
		{
			
			return vector3d<T>(this->data[1] * Vector.data[2] - this->data[2] * Vector.data[1], this->data[2] * Vector.data[0] - this->data[0] * Vector.data[2], this->data[0] * Vector.data[1] - this->data[1] * Vector.data[0]);
		}


		vector3d<T> normalize()
		{
			float Mag = this->Magnitude();
			return vector3d(this->data[0] / Mag, this->data[1]/ Mag, this->data[2] / Mag);
		}



	};

	typedef vector3d<int> Vector3D_int;
	typedef vector3d<float> Vector3D_flo;
	typedef vector3d<double> Vector3D_dou;
	
	template <typename T>
	struct vector4d : public s_matrix<T,4,1>
	{

		vector4d()
		{
			this->data[0] = 0;
			this->data[1] = 0;
			this->data[2] = 0;
			this->data[3] = 0;
		}

		vector4d(T _X, T _Y, T _Z, T _W)
		{

			this->data[0] = _X;
			this->data[1] = _Y;
			this->data[2] = _Z;
			this->data[3] = _W;
		}

		vector4d(vector3d<T> vec)
		{
			this->data[0] = vec.data[0];
			this->data[1] = vec.data[1];
			this->data[2] = vec.data[2];
			this->data[3] = 1;
		}


	};


	typedef vector4d<float> Vector4D_flo;
	typedef vector4d<double> Vector4D_dou;

	struct _Quaternion
	{
		float scalar;
		Vector3D_flo vector;

		_Quaternion()
		{
			scalar = 0;
		}

		_Quaternion(float s, Vector3D_flo vec)
		{
			scalar = s;
			vector = vec;
		}

		_Quaternion(float w, float x, float y, float z)
		{
			scalar = w;
			vector = Vector3D_flo(x, y, z);
		}


		float norm()
		{
			return 1 / sqrtf(scalar * scalar + vector.Dot_product(vector));

		}

		_Quaternion normalize()
		{
			float n = norm();
			return _Quaternion(scalar * n, vector * n);
		}


		_Quaternion operator+(_Quaternion q2)
		{
			return _Quaternion(scalar + q2.scalar, vector + q2.vector);

		}

		_Quaternion operator*(float m)
		{
			return _Quaternion(scalar * m, vector * m);

		}


		_Quaternion operator*(_Quaternion q)
		{
			return _Quaternion(scalar * q.scalar - vector.Dot_product(q.vector), (q.vector * scalar) + (vector * q.scalar) + vector.cross_product(q.vector));

		}


		_Quaternion operator-(_Quaternion q2)
		{
			return _Quaternion(scalar - q2.scalar, vector - q2.vector);

		}

		_Quaternion conjugate()
		{
			return _Quaternion(scalar, -vector);
		}

		_Quaternion inverse()
		{
			float n = norm();

			return conjugate()*n*n;
		}

		friend
			std::ostream& operator<<(std::ostream& os, _Quaternion q)
		{

			os << "[" << q.scalar << "," << q.vector.data[0] << "," << q.vector.data[1] << "," << q.vector.data[2] << "]\n";

			return os;
		}
		

	};
	
	

	_Quaternion axisAngle_quaternion(float angle, float x, float y, float z)
	{

		float q0 = cosf(angle / 2);
		float sinus = sinf(angle / 2);
		float q1 = sinus * x;
		float q2 = sinus * y;
		float q3 = sinus * z;
		return { q0,q1,q2,q3 };

	}
	

	template <int _size>
	struct shape	
	{
		s_matrix<float, 3, _size> positions;
		
		
		s_matrix<float,3,1> get_center()
		{
			s_matrix<float, _size, 1> t;
			for (int i = 0; i < _size; i++)
				t.data[i] = 1;
			
			return ( positions*t)/_size;
		}

		void translate(Vector3D_flo v3)
		{
			for(int i=0;i<3;i++)
				for (int j = 0; j < _size; j++)
				{
					positions.data[i * _size + j] += v3.data[i];
				}

		}

		void rotateX(float angle)
		{
			m_3xRot X(angle);
			positions = X * positions;
		} 

		void rotateZ(float angle)
		{
			m_3zRot Z(angle);
			positions = Z * positions;
		}

		void rotateY(float angle)
		{
			m_3yRot Y(angle);
			positions = Y*positions;
		}

		void rotate(s_matrix<float, 3, 3> mat)
		{
			positions = mat * positions;
		}
		


		
		

	};

	typedef shape<3> triangle;
	

