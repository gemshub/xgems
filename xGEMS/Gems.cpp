// xGEMS is a C++ and Python library for thermodynamic modeling by Gibbs energy minimization
//
// Copyright (C) 2018 Allan Leal, Dmitrii Kulik
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "Gems.hpp"

// C++ includes
#include <chrono>
#include <map>
#include <set>
#include <vector>

// GEMS3K includes
#define IPMGEMPLUGIN
#define NOPARTICLEARRAY
#include <xGEMS/GEMS3K/node.h>

// xGEMS includes
//#include <Reaktoro/Common/Constants.hpp>
//#include <Reaktoro/Common/Exception.hpp>
//#include <Reaktoro/Common/TimeUtils.hpp>
//#include <Reaktoro/Core/ChemicalState.hpp>
//#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace xGEMS {

struct Gems::Impl
{
    /// The TNode instance from Gems
    std::shared_ptr<TNode> node;

    /// The elapsed time of the equilibrate method (in units of s)
    double elapsed_time = 0;

    /// The options for Gems
    GemsOptions options;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a default Impl instance
    Impl(std::string filename)
    {
        // Initialize the GEMS `node` member
        node = std::make_shared<TNode>();
        if(node->GEM_init(filename.c_str()))
            throw std::runtime_error("\n*** ERROR ***\n"
                "Could not initialize the ChemicalEngine object.\n"
                    "Make sure the provided file exists relative to the working directory.");
    }
};

Gems::Gems()
: pimpl(new Impl())
{
}

Gems::Gems(std::string filename)
: pimpl(new Impl(filename))
{
}

Gems::~Gems()
{}

auto Gems::temperature() const -> double
{
    return node()->Get_TK();
}

auto Gems::pressure() const -> double
{
    return node()->Get_P();
}

auto Gems::elementAmounts() const -> VectorConstRef
{
    return Vector::Map(node()->pCNode()->bIC, numElements());
}

auto Gems::speciesAmounts() const -> VectorConstRef
{
    return Vector::Map(node()->pCNode()->xDC, numSpecies());
}

auto Gems::numElements() const -> unsigned
{
    return node()->pCSD()->nIC;
}

auto Gems::numSpecies() const -> unsigned
{
    return node()->pCSD()->nDC;
}

auto Gems::numPhases() const -> unsigned
{
    return node()->pCSD()->nPH;
}

auto Gems::numSpeciesInPhase(Index iphase) const -> unsigned
{
    return node()->pCSD()->nDCinPH[iphase];
}

auto Gems::elementName(Index ielement) const -> std::string
{
    return node()->pCSD()->ICNL[ielement];
}

auto Gems::elementMolarMass(Index ielement) const -> double
{
    return node()->ICmm(ielement);
}

auto Gems::elementStoichiometry(Index ispecies, Index ielement) const -> double
{
    return node()->DCaJI(ispecies, ielement);
}

auto Gems::speciesName(Index ispecies) const -> std::string
{
    return node()->pCSD()->DCNL[ispecies];
}

auto Gems::phaseName(Index iphase) const -> std::string
{
    return node()->pCSD()->PHNL[iphase];
}

auto Gems::setOptions(const GemsOptions& options) -> void
{
    pimpl->options = options;
}

auto Gems::equilibrate(double T, double P, VectorConstRef b) -> void
{
    // Begin timing
    auto begin = std::chrono::high_resolution_clock::now();

    // Set temperature and pressure
    node()->setTemperature(T);
    node()->setPressure(P);

    // Set the mole amounts of the elements
    for(unsigned i = 0; i < numElements(); ++i)
        node()->pCNode()->bIC[i] = b[i];

    // Solve the equilibrium problem with gems
    node()->pCNode()->NodeStatusCH =
        pimpl->options.warmstart ? NEED_GEM_SIA : NEED_GEM_AIA;
    node()->GEM_run(false);

    // Finish timing
    auto end = std::chrono::high_resolution_clock::now();

    // Set the elapsed time member
    pimpl->elapsed_time = std::chrono::duration<double>(end - begin).count();
}

auto Gems::converged() const -> bool
{
    const auto status = node()->pCNode()->NodeStatusCH;
    return status == OK_GEM_AIA || status == OK_GEM_SIA;
}

auto Gems::numIterations() const -> unsigned
{
    return node()->pCNode()->IterDone;
}

auto Gems::elapsedTime() const -> double
{
    return pimpl->elapsed_time;
}

auto Gems::node() const -> std::shared_ptr<TNode>
{
    return pimpl->node;
}

} // namespace xGEMS
