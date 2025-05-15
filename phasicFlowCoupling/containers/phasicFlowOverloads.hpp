#ifndef __phasicFlowOverloads_hpp__
#define __phasicFlowOverloads_hpp__

#include <vector>

#include "iOstream.hpp"
#include "token.hpp"

namespace pFlow::coupling
{

	template<typename T>
	iOstream& operator<< (iOstream& os, const std::vector<T> vec)
	{

	    token::punctuationToken sep;

	    if(vec.size()<10 )
	        sep = token::SPACE;
	    else
	        sep = token::NL;

	    os<<token::BEGIN_LIST;
	    for(const auto& v:vec)
	    {
	        os<<v<<sep;
	    }
	    os<<token::END_LIST;
	    return os;
	}

}


#endif //__phasicFlowOverloads_hpp__


