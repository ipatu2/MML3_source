#pragma once
#include<memory>
#include<cstring>
//#include<type_traits>
#include<initializer_list>


namespace MML3{

	template<typename T>
	class array_
	{
		//static_assert(std::is_trivially_copy_constructible<T>::value, "array_ is suited only for arithmetic types");
	public:
		// CTOR's  ----------------------
		array_() = default;
		array_(const array_& rhs)				{ copy_from(rhs.p_, rhs.sz_); }
		array_(array_&& rhs)					{ swap(rhs); }
		explicit array_(size_t sz)				:p_(new T[sz]), sz_(sz){}
		template<typename S>
		array_(const S* src, size_t sz)			{ copy_from(src,sz); }
		template<typename S>
		array_(std::initializer_list<S> s)		{ copy_from(s.begin(), s.size()); }

		~array_() { free(); }
		// METHODs ----------------------
		void		resize(size_t sz)			{ if (sz_ != sz){ array_ tmp(sz); swap(tmp);} }
		void		free()						{ delete[]p_; p_ = nullptr; sz_ = 0; }
		array_&		swap(array_& rhs)			{ std::swap(p_, rhs.p_); std::swap(sz_, rhs.sz_); return *this; }
		size_t		size()const					{ return sz_; }

		T*			begin()						{return sz_?p_:nullptr;}
		const T*	begin()const				{return sz_?p_:nullptr;}
		T*			end()						{return sz_?p_+sz_:nullptr;}
		const T*	end()const					{return sz_?p_+sz_:nullptr;}

		T&			operator[](size_t i)		{ return p_[i]; }
		const T&	operator[](size_t i)const	{ return p_[i]; }
		array_&     operator=(const array_& rhs){ return copy_from(rhs.p_, rhs.sz_); }
		array_&     operator=(array_&& rhs)		{ free(); std::swap(p_, rhs.p_); std::swap(sz_, rhs.sz_);return *this; }
		template<typename S>
		array_&     operator=(std::initializer_list<S> s)    { return copy_from(s.begin(), s.size()); }
		array_&		fill(T val)					{ for (size_t i = 0; i != sz_; ++i) p_[i]=val; return *this; }

		array_&     copy_from(const T* src, size_t sz )
		{
			resize(sz);
			std::memcpy( (void*)p_, (const void*)src, sizeof(T)*sz);
			return *this;
		}

	private:
		T		*p_  =nullptr;
		size_t	sz_ = 0;
	};


}// end namespace
