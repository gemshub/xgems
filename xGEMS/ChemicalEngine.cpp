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

#include "ChemicalEngine.hpp"

// C++ includes
#include <chrono>

// GEMS3K includes
#define IPMGEMPLUGIN
#define NOPARTICLEARRAY
#include <xGEMS/GEMS3K/node.h>

namespace xGEMS {

struct ChemicalEngine::Impl
{
    /// The TNode instance from ChemicalEngine
    std::shared_ptr<TNode> node;

    /// The elapsed time of the equilibrate method (in units of s)
    double elapsed_time = 0;

    /// The options for ChemicalEngine
    ChemicalEngineOptions options;

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

ChemicalEngine::ChemicalEngine()
: pimpl(new Impl())
{
}

ChemicalEngine::ChemicalEngine(const ChemicalEngine& other)
: pimpl(new Impl(*other.pimpl))
{
}

ChemicalEngine::ChemicalEngine(std::string filename)
: pimpl(new Impl(filename))
{
}

ChemicalEngine::~ChemicalEngine()
{}

auto ChemicalEngine::operator=(ChemicalEngine other) -> ChemicalEngine&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ChemicalEngine::temperature() const -> double
{
    return pimpl->node->Get_TK();
}

auto ChemicalEngine::pressure() const -> double
{
    return pimpl->node->Get_P();
}

auto ChemicalEngine::elementAmounts() const -> VectorConstRef
{
    return Vector::Map(pimpl->node->pCNode()->bIC, numElements());
}

auto ChemicalEngine::speciesAmounts() const -> VectorConstRef
{
    return Vector::Map(pimpl->node->pCNode()->xDC, numSpecies());
}

auto ChemicalEngine::numElements() const -> unsigned
{
    return pimpl->node->pCSD()->nIC;
}

auto ChemicalEngine::numSpecies() const -> unsigned
{
    return pimpl->node->pCSD()->nDC;
}

auto ChemicalEngine::numPhases() const -> unsigned
{
    return pimpl->node->pCSD()->nPH;
}

auto ChemicalEngine::numSpeciesInPhase(Index iphase) const -> unsigned
{
    return pimpl->node->pCSD()->nDCinPH[iphase];
}

auto ChemicalEngine::elementName(Index ielement) const -> std::string
{
    return pimpl->node->pCSD()->ICNL[ielement];
}

auto ChemicalEngine::elementMolarMass(Index ielement) const -> double
{
    return pimpl->node->ICmm(ielement);
}

auto ChemicalEngine::elementStoichiometry(Index ispecies, Index ielement) const -> double
{
    return pimpl->node->DCaJI(ispecies, ielement);
}

auto ChemicalEngine::speciesName(Index ispecies) const -> std::string
{
    return pimpl->node->pCSD()->DCNL[ispecies];
}

auto ChemicalEngine::phaseName(Index iphase) const -> std::string
{
    return pimpl->node->pCSD()->PHNL[iphase];
}

auto ChemicalEngine::setOptions(const ChemicalEngineOptions& options) -> void
{
    pimpl->options = options;
}

auto ChemicalEngine::equilibrate(double T, double P, VectorConstRef b) -> void
{
    // Begin timing
    auto begin = std::chrono::high_resolution_clock::now();

    // Set temperature and pressure
    pimpl->node->setTemperature(T);
    pimpl->node->setPressure(P);

    // Set the mole amounts of the elements
    for(unsigned i = 0; i < numElements(); ++i)
        pimpl->node->pCNode()->bIC[i] = b[i];

    // Solve the equilibrium problem with gems
    pimpl->node->pCNode()->NodeStatusCH =
        pimpl->options.warmstart ? NEED_GEM_SIA : NEED_GEM_AIA;
    pimpl->node->GEM_run(false);

    // Finish timing
    auto end = std::chrono::high_resolution_clock::now();

    // Set the elapsed time member
    pimpl->elapsed_time = std::chrono::duration<double>(end - begin).count();
}

auto ChemicalEngine::converged() const -> bool
{
    const auto status = pimpl->node->pCNode()->NodeStatusCH;
    return status == OK_GEM_AIA || status == OK_GEM_SIA;
}

auto ChemicalEngine::numIterations() const -> unsigned
{
    return pimpl->node->pCNode()->IterDone;
}

auto ChemicalEngine::elapsedTime() const -> double
{
    return pimpl->elapsed_time;
}

} // namespace xGEMS
