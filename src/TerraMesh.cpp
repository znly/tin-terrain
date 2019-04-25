#include "tntn/TerraMesh.h"
#include "tntn/TerraUtils.h"
#include "tntn/logging.h"
#include "tntn/SurfacePoints.h"
#include "tntn/DelaunayTriangle.h"

#include <iostream>
#include <fstream>
#include <array>
#include <unordered_map>
#include <cmath>

namespace tntn {
namespace terra {

void TerraMesh::greedy_insert(double max_error)
{
    m_max_error = max_error;
    m_counter = 0;
    int w = m_raster->get_width();
    int h = m_raster->get_height();
    TNTN_ASSERT(w > 0);
    TNTN_ASSERT(h > 0);

    TNTN_LOG_INFO("starting greedy insertion with raster width: {}, height: {}", w, h);

    // Initialize m_used
    m_used.allocate(w, h);
    m_used.set_all(0);

    // Ensure the four corners are not NAN, otherwise the algorithm can't proceed.
    this->repair_point(0, 0);
    this->repair_point(0, h - 1);
    this->repair_point(w - 1, h - 1);
    this->repair_point(w - 1, 0);

    // Initialize the mesh to two triangles with the height field grid corners as vertices
    TNTN_LOG_INFO("initialize the mesh with four corner points");
    this->init_mesh(
        glm::dvec2(0, 0), glm::dvec2(0, h - 1), glm::dvec2(w - 1, h - 1), glm::dvec2(w - 1, 0));

    m_used.value(0, 0) = 1;
    m_used.value(h - 1, 0) = 1;
    m_used.value(h - 1, w - 1) = 1;
    m_used.value(0, w - 1) = 1;

    // Initialize m_token
    m_token.allocate(w, h);
    m_token.set_all(0);

    // Scan all the triangles and push all candidates into a stack
    dt_ptr t = m_first_face;
    while(t)
    {
        scan_triangle(t);
        t = t->getLink();
    }

    // Iterate until the error threshold is met
    while(!m_candidates.empty())
    {
        Candidate candidate = m_candidates.grab_greatest();

        if(candidate.importance < m_max_error * (candidate.edge ? 0.5 : 1.0)) continue;

        // Skip if the candidate is not the latest
        if(m_token.value(candidate.y, candidate.x) != candidate.token) continue;

        m_used.value(candidate.y, candidate.x) = 1;

        //TNTN_LOG_DEBUG("inserting point: ({}, {}, {})", candidate.x, candidate.y, candidate.z);
        this->insert(glm::dvec2(candidate.x, candidate.y), candidate.triangle);
    }

    TNTN_LOG_INFO("finished greedy insertion");
}

void TerraMesh::scan_triangle_line(const Plane& plane,
                                   int y,
                                   double x1,
                                   double x2,
                                   Candidate& candidate,
                                   const double no_data_value)
{
    const int startx = static_cast<int>(ceil(fmin(x1, x2)));
    const int endx = static_cast<int>(floor(fmax(x1, x2)));

    if(startx > endx) return;

    double z0 = plane.eval(startx, y);
    double dz = plane.a;

    for(int x = startx; x <= endx; x++)
    {
        if(!m_used.value(y, x))
        {
            const double z = m_raster->value(y, x);
            if(!is_no_data(z, no_data_value))
            {
                const double diff = fabs(z - z0);
                //TNTN_LOG_DEBUG("candidate consider: ({}, {}, {}), diff: {}", x, y, z, diff);
                candidate.consider(x, y, z, diff);
            }
        }
        z0 += dz;
    }
}

void TerraMesh::scan_line_vertical(const Plane& plane,
                                   int x,
                                   double y1,
                                   double y2,
                                   Candidate& candidate,
                                   const double no_data_value)
{
    const int starty = static_cast<int>(ceil(fmin(y1, y2)));
    const int endy = static_cast<int>(floor(fmax(y1, y2)));

    if(starty > endy) return;

    double z0 = plane.eval(x, starty);
    double dz = plane.b;

    for(int y = starty; y <= endy; y++)
    {
        if(!m_used.value(y, x))
        {
            const double z = m_raster->value(y, x);
            if(!is_no_data(z, no_data_value))
            {
                const double diff = fabs(z - z0);
                //TNTN_LOG_DEBUG("candidate consider: ({}, {}, {}), diff: {}", x, y, z, diff);
                candidate.consider(x, y, z, diff);
            }
        }
        z0 += dz;
    }
}

void TerraMesh::scan_triangle(dt_ptr t)
{
    Plane z_plane;
    compute_plane(z_plane, t, *m_raster);

    std::array<Point2D, 3> by_y = {{
        t->point1(),
        t->point2(),
        t->point3(),
    }};
    order_triangle_points(by_y);
    const double v0_x = by_y[0].x;
    const double v0_y = by_y[0].y;
    const double v1_x = by_y[1].x;
    const double v1_y = by_y[1].y;
    const double v2_x = by_y[2].x;
    const double v2_y = by_y[2].y;

    const double no_data_value = m_raster->get_no_data_value();

    Candidate candidate_v = {0, 0, 0.0, -DBL_MAX, m_counter++, t};
    Candidate candidate_h = {0, 0, 0.0, -DBL_MAX, m_counter++, t};
    const double w = m_raster->get_width()-1;
    const double h = m_raster->get_height()-1;

    for (int i = 0; i < 3; i++) {
        int ii = (i + 1) % 3;

        if ((by_y[i].x == 0 && by_y[ii].x == 0) ||
            (by_y[i].x == w && by_y[ii].x == w)) {
            //printf("scane v_edge\n");

            // edge left/right
            scan_line_vertical(z_plane, by_y[i].x, by_y[i].y, by_y[ii].y, candidate_v, no_data_value);

        } else if ((by_y[i].y == 0 && by_y[ii].y == 0) ||
                   (by_y[i].y == h && by_y[ii].y == h)) {
            //printf("scan h_edge\n");
            // edge top/bottom
            scan_triangle_line(z_plane, by_y[i].y, by_y[i].x, by_y[ii].x, candidate_h, no_data_value);
        }
    }

    bool v_edge = candidate_v.importance != -DBL_MAX;
    bool h_edge = candidate_h.importance != -DBL_MAX;

    if (v_edge || h_edge) {
        if (v_edge) {
            //printf("v_edge\n");
            m_token.value(candidate_v.y, candidate_v.x) = candidate_v.token;
            candidate_v.edge = true;
            m_candidates.push_back(candidate_v);
        }
        if (h_edge) {
            //printf("h_edge\n");
            m_token.value(candidate_h.y, candidate_h.x) = candidate_h.token;
            candidate_h.edge = true;
            m_candidates.push_back(candidate_h);
        }
        return;
    }

    Candidate candidate = {0, 0, 0.0, -DBL_MAX, m_counter++, t};

    const double dx2 = (v2_x - v0_x) / (v2_y - v0_y);

    if(v1_y != v0_y)
    {
        const double dx1 = (v1_x - v0_x) / (v1_y - v0_y);

        double x1 = v0_x;
        double x2 = v0_x;

        const int starty = v0_y;
        const int endy = v1_y;

        for(int y = starty; y < endy; y++)
        {
            scan_triangle_line(z_plane, y, x1, x2, candidate, no_data_value);
            x1 += dx1;
            x2 += dx2;
        }
    }

    if(v2_y != v1_y)
    {
        const double dx1 = (v2_x - v1_x) / (v2_y - v1_y);

        double x1 = v1_x;
        double x2 = v0_x;

        const int starty = v1_y;
        const int endy = v2_y;

        for(int y = starty; y <= endy; y++)
        {
            scan_triangle_line(z_plane, y, x1, x2, candidate, no_data_value);
            x1 += dx1;
            x2 += dx2;
        }
    }

    // We have now found the appropriate candidate point
    m_token.value(candidate.y, candidate.x) = candidate.token;

    // Push the candidate into the stack
    m_candidates.push_back(candidate);
}

std::unique_ptr<Mesh> TerraMesh::convert_to_mesh()
{
    // Find all the vertices
    int w = m_raster->get_width();
    int h = m_raster->get_height();

    std::vector<Vertex> mvertices;

    Raster<int> vertex_id;
    vertex_id.allocate(w, h);
    vertex_id.set_all(0);

    const double no_data_value = m_raster->get_no_data_value();
    int index = 0;
    for(int y = 0; y < h; y++)
    {
        for(int x = 0; x < w; x++)
        {
            if(m_used.value(y, x) == 1)
            {
                const double z = m_raster->value(y, x);
                if(is_no_data(z, no_data_value))
                {
                    continue;
                }

                Vertex v = Vertex({m_raster->col2x(x), m_raster->row2y(y), z});
                mvertices.push_back(v);
                vertex_id.value(y, x) = index;
                index++;
            }
        }
    }

    // Find all the faces
    std::vector<Face> mfaces;
    dt_ptr t = m_first_face;
    while(t)
    {
        Face f;

        glm::dvec2 p1 = t->point1();
        glm::dvec2 p2 = t->point2();
        glm::dvec2 p3 = t->point3();

        if(!ccw(p1, p2, p3))
        {
            f[0] = vertex_id.value((int)p1.y, (int)p1.x);
            f[1] = vertex_id.value((int)p2.y, (int)p2.x);
            f[2] = vertex_id.value((int)p3.y, (int)p3.x);
        }
        else
        {
            f[0] = vertex_id.value((int)p3.y, (int)p3.x);
            f[1] = vertex_id.value((int)p2.y, (int)p2.x);
            f[2] = vertex_id.value((int)p1.y, (int)p1.x);
        }

        mfaces.push_back(f);

        t = t->getLink();
    }

    // now initialise our mesh class with this
    auto mesh = std::make_unique<Mesh>();
    mesh->from_decomposed(std::move(mvertices), std::move(mfaces));
    return mesh;
}

} //namespace terra
} //namespace tntn
