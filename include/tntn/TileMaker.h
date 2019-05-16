#pragma once

#include "tntn/Mesh.h"
#include "tntn/MeshWriter.h"
#include "tntn/FileFormat.h"

#include <memory>

namespace tntn {

bool write_mesh_to_file(const char* filename, const Mesh& m, const FileFormat& f);

class TileMaker
{
    std::unique_ptr<Mesh> m_mesh;

  public:
    TileMaker() : m_mesh(std::make_unique<Mesh>()) {}

    void setMeshWriter(MeshWriter* w);
    bool loadObj(const char* filename);
    void loadMesh(std::unique_ptr<Mesh> mesh);
    // void dumpTile(int tx, int ty, int zoom, const char* filename);
    bool dumpTile(int tx, int ty, int zoom, const char* filename, MeshWriter& mw);
};

} //namespace tntn
