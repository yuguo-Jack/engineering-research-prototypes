#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace m6 {
struct ScalarSourceLinearization{double diagonal,source;};
// Positive-diagonal convention. Fold a negative deferred source into the
// diagonal while preserving diagonal*previous-source at this linearization.
inline ScalarSourceLinearization linearize_negative_scalar_source(double previous,double diagonal,double source){
 if(!(previous>0&&std::isfinite(previous)&&std::isfinite(diagonal)&&std::isfinite(source)))throw std::runtime_error("invalid scalar source linearization inputs");
 if(source<0){diagonal-=source/previous;source=0.;}
 return {diagonal,source};
}
// Physical scalar units. Clip the absolute solution before blending, and
// retain the computed (1-relaxation), including its ordinary floating rounding.
// No nextafter and no second clipping of the blended result are applied.
inline double clipped_relaxed_scalar(double previous,double absolute_solution,double minimum,double relaxation){
 if(!(previous>0&&minimum>0&&relaxation>0&&relaxation<=1&&std::isfinite(previous)&&std::isfinite(absolute_solution)&&std::isfinite(minimum)))throw std::runtime_error("invalid scalar update inputs");
 return (1-relaxation)*previous+relaxation*std::max(absolute_solution,minimum);
}
}
