/*------------------------------- phasicFlow ---------------------------------
      O        C enter of
     O O       E ngineering and
    O   O      M ultiscale modeling of
   OOOOOOO     F luid flow       
------------------------------------------------------------------------------
  Copyright (C): www.cemf.ir
  email: hamid.r.norouzi AT gmail.com
------------------------------------------------------------------------------  
Licence:
  This file is part of phasicFlow code. It is a free software for simulating 
  granular and multiphase flows. You can redistribute it and/or modify it under
  the terms of GNU General Public License v3 or any other later versions. 
 
  phasicFlow is distributed to help others in their research in the field of 
  granular and multiphase flows, but WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

-----------------------------------------------------------------------------*/

#include "sphereDrag.hpp"
#include "self.hpp"
#include "Gaussian.hpp"
#include "GaussianIntegral.hpp"
#include "adaptiveGaussian.hpp"
#include "diffusion.hpp"

#include "DiFelice.hpp"
#include "Rong.hpp"
#include "ErgunWenYu.hpp"

#define makeSphereDrag(closure)                 \
                                                \
template class pFlow::coupling::sphereDrag      \
    <                                           \
        pFlow::coupling::self,                  \
        closure,                                \
        true                                    \
    >;                                          \
template class pFlow::coupling::sphereDrag      \
    <                                           \
        pFlow::coupling::self,                  \
        closure,                                \
        false                                   \
    >;                                          \
template class pFlow::coupling::sphereDrag      \
    <                                       \
        pFlow::coupling::Gaussian,          \
        closure,                            \
        true                                \
    >;                                      \
template class pFlow::coupling::sphereDrag  \
    <                                       \
        pFlow::coupling::Gaussian,          \
        closure,                            \
        false                               \
    >;                                      \
template class pFlow::coupling::sphereDrag  \
    <                                       \
        pFlow::coupling::GaussianIntegral,  \
        closure,                            \
        true                                \
    >;                                      \
template class pFlow::coupling::sphereDrag  \
    <                                       \
        pFlow::coupling::GaussianIntegral,  \
        closure,                            \
        false                               \
    >;                                      \
template class pFlow::coupling::sphereDrag  \
    <                                       \
        pFlow::coupling::adaptiveGaussian,  \
        closure,                            \
        true                                \
    >;                                      \
template class pFlow::coupling::sphereDrag  \
    <                                       \
        pFlow::coupling::adaptiveGaussian,  \
        closure,                            \
        false                               \
    >;                                      \
template class pFlow::coupling::sphereDrag  \
    <                                       \
        pFlow::coupling::diffusion,         \
        closure,                            \
        true                                \
    >;                                      \
template class pFlow::coupling::sphereDrag  \
    <                                       \
        pFlow::coupling::diffusion,         \
        closure,                            \
        false                               \
    >;                                                         

makeSphereDrag(pFlow::coupling::DiFelice);
makeSphereDrag(pFlow::coupling::ErgunWenYu);
makeSphereDrag(pFlow::coupling::Rong);
