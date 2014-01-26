
#pragma once
#include"MML3-config.h"
#include"MML3-Types.h"
#include"MML3-Matrix.h"
#include<valarray>
#include<memory>
#if defined(MML3_LAPACK)
	#include"MML3-lapack.h"
#else
#error MML3::Solver relies upon Lapacke
#endif



namespace MML3{

	namespace Solver{


		// errors codes in the range (-9, -1) are reserved for lapack errors : -i means that the i'th parameter passed to lapack routine is wrong
		// error codes in the range (1,...)  in the factorization routine indicates that the i'th equation cannot be factorized
		// 0 return codes means that its all right
		enum State { OK = 0, NO_INIT = -100, BAD_SZ = -10, BAD_ON_SOLVE = -11, BAD_MTYPE };

		template<typename MAT>
		class LU
		{
		public:
			typedef typename MAT::matrix_t		matrix_t;
			typedef typename MAT::value_t 		value_t;
			typedef	lapack_int_t				int_t;
			static_assert(MAT::GE && MAT::RE, " only General Rectangular matrices are supported by LU solver");
			
			
		public:
			LU(matrix_t& A){ factorize(A); }
			// factorizes the matrix A and mantains internally a reference to A and a pivot vector
			~LU(){pA_ = nullptr;}


			bool	factorize(matrix_t& A);
			bool	solve(matrix_t& X, bool trA = false);
			int		state()const{return state_;}
			bool	good()const{ return state_ == State::OK; }

			operator bool()const{ return good();}
		private:
			int state_ = NO_INIT;
			std::valarray<int_t> pivot_;
			matrix_t*   pA_ = nullptr;
		};

		

		//  LLT solver for symmetric positive definite  matrices 
		template<typename MAT>
		class LLT
		{
			typedef typename MAT::matrix_t		matrix_t;
			typedef typename MAT::value_t 		value_t;
			typedef	lapack_int_t				int_t;


			//static_assert(MAT::IS::RE || (MAT::IS::Packed  && MAT::IS::ColMajor),
			//	"Packed Symmetric Row major matrices are not fully supported by Lapacke ?pptrf\s function (version 11.1)");
			// but i provided a workaround
			
		public:
			LLT(matrix_t& A, bool L = true) { factorize(A, L); }
			~LLT(){ pA_ = nullptr; }

			int  state()const{ return state_; }
			bool good()const{ return state_ == 0;}
			operator bool()const{ return good(); }

			bool factorize(matrix_t& A, bool L = true);

			template<typename BMAT>
			bool solve(BMAT& X);
		private:
			int state_ = NO_INIT;
			matrix_t*  pA_=nullptr;
			char uplo_;

		};

		

		


		
		//////////////////////////////////////////////
		// IMPLEMENTATION


		template<typename MAT>
		bool LU<MAT>::factorize(matrix_t& A)
		{
			pA_ = &A;
			pivot_.resize(A.nrows());
			state_=LAPACK<value_t>::getrf(A.CblasOrder(), A.nrows(), A.ncols(), A.begin(), A.leading_dim(), &pivot_[0]);
			return good();
		}




		template<typename MAT>
		bool LU<MAT>::solve(matrix_t& X, bool trA = false)
		{
			if (state_ != 0)
				state_ = State::BAD_ON_SOLVE;
			else if (X.nrows() != pA_->nrows())
				state_ = State::BAD_SZ;
			else
			{
				char tr('N');
				if (trA)
					tr = 'S';
				state_ = LAPACK<value_t>::getrs(pA_->CblasOrder(), tr, pA_->nrows(), X.ncols(), pA_->begin(), pA_->leading_dim(), &pivot_[0], X.begin(), X.leading_dim());
			}
			return good();
		}




		/*template<typename MAT>
		template<>
		bool LLT<MAT>::factorize_<MAT_SHAPE::RE, MAT_ORD::COL>(matrix_t& A, bool L )
		{
			char uplo_ = L ? 'L' : 'U';
			pA_ = &A;
			
			state_ = LAPACK<value_t >::potrf(A.CblasOrder(), uplo_, A.nrows(), A.begin(), A.leading_dim());
			return good();
		}

		template<typename MAT>
		template<>
		bool LLT<MAT>::factorize_<MAT_SHAPE::RE, MAT_ORD::ROW>(matrix_t& A, bool L)
		{
			char uplo_ = L ? 'L' : 'U';
			pA_ = &A;

			state_ = LAPACK<value_t >::potrf(A.CblasOrder(), uplo_, A.nrows(), A.begin(), A.leading_dim());
			return good();
		}*/




		template<typename MAT>
		bool LLT<MAT>::factorize(matrix_t& A, bool L )
		{
			uplo_ = L ? 'L' : 'U';
			pA_ = &A;
			// for rectangular matrices there is no problem for lapacke potrf
			if (A.IS::RE)
				state_=LAPACK<value_t >::potrf(A.CblasOrder(), uplo_, A.nrows(), A.begin(), A.leading_dim());
			// while for packed matrices it seems that intel mkl pptrf does not handle properly row major formats, so ...
			else
			{
				// Col major formats are treated correctly
				if (A.IS::ColMajor)
				{
					uplo_ = A.IS::LT?'L':'U';
					state_ = LAPACK<value_t >::pptrf(A.CblasOrder(), uplo_, A.nrows(), A.begin());
				}
				// Row major formats require some workaround
				else
				{
					// I treat Row major  packed LT mats as they where ColMajor packed UT because the underlying linear vector is the same
					if (A.IS::LT)
					{
						uplo_ = 'L';
						state_ = LAPACK<value_t >::pptrf(Option::ColMajor, 'U', A.nrows(), A.begin());
					}
					// Row major packed UT have the same array structures of Col major LT
					else
					{
						uplo_ = 'U';
						state_ = LAPACK<value_t >::pptrf(Option::ColMajor, 'L', A.nrows(), A.begin());
					}
				}
			}
			return good();
		}

		template<typename MAT>
		template<typename BMAT>
		bool LLT<MAT>::solve(BMAT& X)
		{
			if (state_ != 0)
				state_= State::BAD_ON_SOLVE;
			else if (X.nrows() != pA_->nrows())
				state_=State::BAD_SZ;
			else
			{
				if (pA_->IS::RE)
					state_ = LAPACK<value_t>::potrs(pA_->CblasOrder(), uplo_, pA_->nrows(), X.ncols(), pA_->begin(), pA_->leading_dim(), X.begin(), X.leading_dim());
				else
				{
					if (pA_->IS::ColMajor)
						state_ = LAPACK<value_t>::pptrs(pA_->CblasOrder(), uplo_, pA_->nrows(), X.ncols(), pA_->begin(), X.begin(), X.leading_dim());
					else
					{
						char tmp_uplo = pA_->IS::LT?'U':'L';
						
						// i make a transposed copy of the rhs
						auto B = transpose(X);
						state_ = LAPACK<value_t>::pptrs(Option::ColMajor, tmp_uplo, pA_->nrows(), X.ncols(), pA_->begin(), B.begin(), B.ncols());
						X = transpose(B);
					}
				}
			}
			return good();
		}
		
		/*

		//  LLT solver specialized for Packed symmetric  matrices 
		template<typename T>
		class LLT<Matrix<T, MAT_PROP::SYM>>
		{
			typedef Matrix<T, MAT_PROP::SYM>			sym_matrix_t;
			typedef typename sym_matrix_t::property_t		property_t;
			typedef typename sym_matrix_t::value_t 			value_t;
			typedef typename LAPACK<value_t>::int_t		int_t;
			typedef Matrix<T, MAT_PROP::GE>				rhs_matrix_t;
		public:
			LLT() = default;
			int factorize(sym_matrix_t& A)
			{
				pA_ = &A;
				// !!! con ROW_MAJOR ed 'L' non funziona e pertanto, in maniera equivalente, chiamo la routine
				// con COL_MAJOR e 'U'
				int info = LAPACK<T>::pptrf(MML3_LAPACK_COL_MAJOR, 'U', A.nrows(), A.begin());
				return info;
			}
			int solve(rhs_matrix_t& X)
			{
				if (X.nrows() != pA_->nrows())
					return std::numeric_limits<int>::min();
				// pare che pptrs non voglia saperne di funzionare in modalità ROW_MAJOR, pertanto eseguo la trasposta del termine noto
				X.set_transpose();
				T* px = X.begin();
				int info = LAPACK<T>::pptrs(MML3_LAPACK_COL_MAJOR, 'U', pA_->nrows(), X.nrows(), pA_->begin(), px, X.leading_dimension());
				X.set_transpose();
				return info;
			}
		private:
			sym_matrix_t*  pA_{};
		};

		// an LLT maker
		template< typename T, typename MP>
		auto  make_LLT(const Matrix<T, MP>&)->LLT< Matrix<T, MP>>
		{
			return LLT<Matrix<T, MP>>();
		};
		*/
		
		

	} // end namespace
} // end namespace
