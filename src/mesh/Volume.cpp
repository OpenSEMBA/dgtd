
#include <metis.h>
#if METIS_VER_MAJOR < 5
#error "Mesh partitioning requires METIS version 5+"
#endif

#include "Volume.h"

namespace SEMBA {
namespace dgtd {
namespace Mesh {

Volume::Volume(const Geometry::Mesh::Unstructured& uns) {
    this->coords() = *uns.coords().clone();
    this->elems() = *uns.elems().clone();
    this->layers() = *uns.layers().clone();
}

vector<vector<Geometry::ElemId>> Volume::getPartitionsIds(
        const size_t nDivisions,
        const vector<pair<Geometry::ElemId,int>> idWgt,
        const Math::Real* taskPower) const {
    // Metis v5 manual:
    // [...] take as input the element-node array of the mesh and
    // compute a k-way partitioning for both its elements and its nodes
    // idWgt contains id and weight pairs.
    vector<vector<Geometry::ElemId>> res;
    res.resize(nDivisions, vector<Geometry::ElemId>());
    // Accounts for the one partition case.
    if (nDivisions == 1) {
        Geometry::ConstElemRGroup physVol = elems();
        physVol.removeMatId(MatId(0));
        const size_t nK = physVol.sizeOf<Geometry::VolR>();
        res[0].resize(nK, Geometry::ElemId(0));
        for (size_t i = 0; i < nK; i++) {
            res[0][i] = (elems())(i)->getId();
        }
        return res;
    }
#ifdef MESH_ALLOW_PARTITIONING
    // Prepares mesh info.
    cout << " - Preparing mesh info... " << flush;
    idx_t ne = elems().sizeOf<Geometry::VolR>();
    idx_t *eptr, *eind;
    eptr = new idx_t[ne+1];
    eind = new idx_t[ne*4];
    size_t counter = 0;
    eptr[0] = counter;
    for (idx_t i = 0; i < ne; i++) {
        const Geometry::VolR* vol = elem_.tet[i];
        for (size_t j = 0; j < vol->numberOfVertices(); j++) {
            eind[counter++] = vol->getVertex(j)->id - 1;
        }
        eptr[i+1] = counter;
    }
    cout << "OK" << endl;
    // Relabels ids, needed by quadratic or linearized meshes.
    cout << " - Relabeling... " << flush;
    DynMatrix<Math::Int> id(ne*4,3);
    for (Math::Int i = 0; i < ne*4; i++) {
        id(i,0) = i;
        id(i,1) = eind[i];
        id(i,2) = 0;
    }
    id.sortRows_omp(1,1);
    Math::Int label = 0;
    for (Math::Int i = 1; i < ne*4; i++) {
        if (id(i,1) == id(i-1,1)) {
            id(i,2) = label;
        } else {
            id(i,2) = ++label;
        }
    }
    id.sortRows_omp(0,0);
    for (Math::Int i = 0; i < ne*4; i++) {
        eind[i] = id(i,2);
    }
    idx_t nn = label+1; // Number of vertices.
    cout << "OK" << endl;
    // Copies weights.
    cout << " - Copying weights... " << flush;
    idx_t *vwgt;
    if (idWgt.size() == 0) {
        vwgt = NULL;
    } else {
        vwgt = new idx_t[ne];
        for (Math::Int e = 0; e < ne; e++) {
            vwgt[e] = idWgt[e].second;
        }
    }
    idx_t *vsize = NULL;
    idx_t nparts = nDivisions;
    idx_t objval;
    idx_t *epart;
    epart = new idx_t[ne];
    idx_t *npart;
    npart = new idx_t[nn];
    cout << "OK" << endl;
    // Computes task computational powers.
    real_t *tpwgts = NULL;
    if (taskPower != NULL) {
        tpwgts = new real_t[nDivisions];
        real_t sum = 0.0;
        for (size_t i = 0; i < nDivisions; i++) {
            tpwgts[i] = taskPower[i];
            sum += tpwgts[i];
        }
        assert(std::abs(sum) - 1.0e-16 < 1.0);
    }
    // METIS options.
    cout << " - Setting Options... " << flush;
    idx_t options[METIS_NOPTIONS];
    Math::Int status;
    status = METIS_SetDefaultOptions(options);
    options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    options[METIS_OPTION_SEED] = (idx_t) 0;
    //	options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
    // c numbering. Starts from 0.
    options[METIS_OPTION_NUMBERING] = 0;
    cout << "OK" << endl;
    // Calls METIS partition function for meshes.
    idx_t ncommon = 3; // Number of common vertices per element.
    cout << " - Calling Part Mesh Dual... " << flush;
    status = METIS_PartMeshDual(
            &ne, &nn, eptr, eind, vwgt, vsize, &ncommon, &nparts,
            tpwgts, options, &objval, epart, npart);
    if (status != METIS_OK) {
        throw Error("METIS_PartMeshDual fn failed with error: " + status);
    }
    cout << "OK" << endl;
    // Converts result.
    for (size_t i = 0; i < nDivisions; i++) {
        res[i].reserve(ne);
    }
    for (Math::Int i = 0; i < ne; i++) {
        size_t id = elem_.tet[i]->getId();
        res[epart[i]].push_back(id);
    }
    // Frees memory.
    delete vwgt;
    delete epart;
    delete npart;
    delete eptr;
    delete eind;
    // Returns result.
    return res;
#else
    throw logic_error("Mesh partitioning is not allowed.");
#endif
}

}
}
}
