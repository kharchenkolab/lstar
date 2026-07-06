// SPDX-License-Identifier: MIT

#ifndef LIBZARR_LIBZARR_HPP
#define LIBZARR_LIBZARR_HPP

/// \file libzarr.hpp
/// Umbrella header: includes the whole libzarr core. Adapters (which may
/// require OS facilities the core must not depend on) are never included
/// here; include them individually from libzarr/adapters/.

#include "libzarr/array.hpp"
#include "libzarr/codecs.hpp"
#include "libzarr/group.hpp"
#include "libzarr/metadata.hpp"
#include "libzarr/sharding.hpp"
#include "libzarr/store.hpp"
#include "libzarr/types.hpp"
#include "libzarr/v2.hpp"
#include "libzarr/v3.hpp"
#include "libzarr/zip.hpp"

#endif  // LIBZARR_LIBZARR_HPP
