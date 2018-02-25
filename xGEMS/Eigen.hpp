// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#pragma once

// Eigen includes
#include <eigen3/Eigen/Core>

namespace xGEMS {

using Vector           = Eigen::VectorXd;                             /// < Alias to Eigen type VectorXd.
using VectorRef        = Eigen::Ref<Eigen::VectorXd>;                 /// < Alias to Eigen type Ref<VectorXd>.
using VectorStridedRef = Eigen::Ref<Vector, 0, Eigen::InnerStride<>>; /// < Alias to Eigen type Ref<VectorXd>.
using VectorConstRef   = Eigen::Ref<const Eigen::VectorXd>;           /// < Alias to Eigen type Ref<const VectorXd>.
using VectorMap        = Eigen::Map<Eigen::VectorXd>;                 /// < Alias to Eigen type Map<VectorXd>.
using VectorConstMap   = Eigen::Map<const Eigen::VectorXd>;           /// < Alias to Eigen type Map<const VectorXd>.

using RowVector           = Eigen::RowVectorXd;                                      /// < Alias to Eigen type RowVectorXd.
using RowVectorRef        = Eigen::Ref<Eigen::RowVectorXd>;                          /// < Alias to Eigen type Ref<RowVectorXd>.
using RowVectorStridedRef = Eigen::Ref<Eigen::RowVectorXd, 0, Eigen::InnerStride<>>; /// < Alias to Eigen type Ref<RowVectorXd>.
using RowVectorConstRef   = Eigen::Ref<const Eigen::RowVectorXd>;                    /// < Alias to Eigen type Ref<const RowVectorXd>.
using RowVectorMap        = Eigen::Map<Eigen::RowVectorXd>;                          /// < Alias to Eigen type Map<VectorXd>.
using RowVectorConstMap   = Eigen::Map<const Eigen::RowVectorXd>;                    /// < Alias to Eigen type Map<const VectorXd>.

using Matrix         = Eigen::MatrixXd;                   ///< Alias to Eigen type MatrixXd.
using MatrixRef      = Eigen::Ref<Eigen::MatrixXd>;       ///< Alias to Eigen type Ref<MatrixXd>.
using MatrixConstRef = Eigen::Ref<const Eigen::MatrixXd>; ///< Alias to Eigen type Ref<const MatrixXd>.
using MatrixMap      = Eigen::Map<Eigen::MatrixXd>;       ///< Alias to Eigen type Map<MatrixXd>.
using MatrixConstMap = Eigen::Map<const Eigen::MatrixXd>; ///< Alias to Eigen type Map<const MatrixXd>.

using Vector = Eigen::VectorXd; /// Alias to Eigen type Eigen::VectorXd.
using VectorXd = Eigen::VectorXd; /// Alias to Eigen type Eigen::VectorXd.
using VectorXi = Eigen::VectorXi; /// Alias to Eigen type Eigen::VectorXi.

using VectorRef = Eigen::Ref<VectorXd>; ///< Alias to Eigen type Eigen::Ref<VectorXd>.
using VectorXdRef = Eigen::Ref<VectorXd>; ///< Alias to Eigen type Eigen::Ref<VectorXd>.
using VectorXiRef = Eigen::Ref<VectorXi>; ///< Alias to Eigen type Eigen::Ref<VectorXi>.

using VectorConstRef = Eigen::Ref<const VectorXd>; ///< Alias to Eigen type Eigen::Ref<const VectorXd>.
using VectorXdConstRef = Eigen::Ref<const VectorXd>; ///< Alias to Eigen type Eigen::Ref<const VectorXd>.
using VectorXiConstRef = Eigen::Ref<const VectorXi>; ///< Alias to Eigen type Eigen::Ref<const VectorXi>.

using VectorMap = Eigen::Map<VectorXd>; ///< Alias to Eigen type Eigen::Map<VectorXd>.
using VectorXdMap = Eigen::Map<VectorXd>; ///< Alias to Eigen type Eigen::Map<VectorXd>.
using VectorXiMap = Eigen::Map<VectorXi>; ///< Alias to Eigen type Eigen::Map<VectorXi>.

using VectorConstMap = Eigen::Map<const VectorXd>; ///< Alias to Eigen type Eigen::Map<const VectorXd>.
using VectorXdConstMap = Eigen::Map<const VectorXd>; ///< Alias to Eigen type Eigen::Map<const VectorXd>.
using VectorXiConstMap = Eigen::Map<const VectorXi>; ///< Alias to Eigen type Eigen::Map<const VectorXi>.

using Matrix = Eigen::MatrixXd; /// Alias to Eigen type Eigen::MatrixXd.
using MatrixXd = Eigen::MatrixXd; /// Alias to Eigen type Eigen::MatrixXd.
using MatrixXi = Eigen::MatrixXi; /// Alias to Eigen type Eigen::MatrixXi.

using MatrixRef = Eigen::Ref<MatrixXd>; ///< Alias to Eigen type Eigen::Ref<MatrixXd>.
using MatrixXdRef = Eigen::Ref<MatrixXd>; ///< Alias to Eigen type Eigen::Ref<MatrixXd>.
using MatrixXiRef = Eigen::Ref<MatrixXi>; ///< Alias to Eigen type Eigen::Ref<MatrixXi>.

using MatrixConstRef = Eigen::Ref<const MatrixXd>; ///< Alias to Eigen type Eigen::Ref<const MatrixXd>.
using MatrixXdConstRef = Eigen::Ref<const MatrixXd>; ///< Alias to Eigen type Eigen::Ref<const MatrixXd>.
using MatrixXiConstRef = Eigen::Ref<const MatrixXi>; ///< Alias to Eigen type Eigen::Ref<const MatrixXi>.

using MatrixMap = Eigen::Map<MatrixXd>; ///< Alias to Eigen type Eigen::Map<MatrixXd>.
using MatrixXdMap = Eigen::Map<MatrixXd>; ///< Alias to Eigen type Eigen::Map<MatrixXd>.
using MatrixXiMap = Eigen::Map<MatrixXi>; ///< Alias to Eigen type Eigen::Map<MatrixXi>.

using MatrixConstMap = Eigen::Map<const MatrixXd>; ///< Alias to Eigen type Eigen::Map<const MatrixXd>.
using MatrixXdConstMap = Eigen::Map<const MatrixXd>; ///< Alias to Eigen type Eigen::Map<const MatrixXd>.
using MatrixXiConstMap = Eigen::Map<const MatrixXi>; ///< Alias to Eigen type Eigen::Map<const MatrixXi>.

} // namespace xGEMS
