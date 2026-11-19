
struct register_embree_ghost_dynamic_op
{
    uint32_t operator()(const Sphere &sphere) const;
    uint32_t operator()(const TriangleMesh &mesh) const;
    const RTCDevice &device;
    const RTCScene &scene;
    const RTCScene &ghost_scene;
    const int time_steps;
    Dynamic &ghost_dynamic;
};

void intersection_filter(const RTCFilterFunctionNArguments *args)
{
    /* avoid crashing when debug visualizations are used */
    if (args->context == nullptr)
        return;
    assert(args->N == 1);
    int *valid = args->valid;
    const IntersectContext *context = (const IntersectContext *)args->context;
    RTCRay *rtc_ray = (RTCRay *)args->ray;
    RTCHit *rtc_hit = (RTCHit *)args->hit;

    /* ignore inactive rays */
    if (valid[0] != -1)
        return;
    valid[0] = 0; //ignore the hit
    //set the payload
    {
        Payload* payload = (Payload*)context->payload;
        payload->dist_t.push_back(rtc_ray->tfar);
        payload->ghost_shape_ids.push_back(rtc_hit->geomID);
        Vector3 position = Vector3{rtc_ray->org_x, rtc_ray->org_y, rtc_ray->org_z} +
                      Vector3{rtc_ray->dir_x, rtc_ray->dir_y, rtc_ray->dir_z} * Real(rtc_ray->tfar);
        Vector3 geometric_normal = normalize(Vector3{rtc_hit->Ng_x, rtc_hit->Ng_y, rtc_hit->Ng_z});
        payload->ghost_hit_positions.push_back(position);
        payload->ghost_hit_geom_normals.push_back(geometric_normal);
        payload->num_ghost_hits += 1;
    }
    
}

void occlusion_filter(const RTCFilterFunctionNArguments *args)
{
    /* avoid crashing when debug visualizations are used */
    if (args->context == nullptr)
        return;

    assert(args->N == 1);
    int *valid = args->valid;
    const IntersectContext *context = (const IntersectContext *)args->context;
    RTCRay *rtc_ray = (RTCRay *)args->ray;
    RTCHit *rtc_hit = (RTCHit *)args->hit;
    /* ignore inactive rays */
    if (valid[0] != -1)
        return;
    valid[0] = 0; //ignore the hit
    //set the payload
    {
        Payload* payload = (Payload*)context->payload;
        payload->dist_t.push_back(rtc_ray->tfar);
        payload->ghost_shape_ids.push_back(rtc_hit->geomID);
        Vector3 position = Vector3{rtc_ray->org_x, rtc_ray->org_y, rtc_ray->org_z} +
                      Vector3{rtc_ray->dir_x, rtc_ray->dir_y, rtc_ray->dir_z} * Real(rtc_ray->tfar);
        Vector3 geometric_normal = normalize(Vector3{rtc_hit->Ng_x, rtc_hit->Ng_y, rtc_hit->Ng_z});
        payload->ghost_hit_positions.push_back(position);
        payload->ghost_hit_geom_normals.push_back(geometric_normal);
        payload->num_ghost_hits += 1;
    }
}
////////////////////////////////

uint32_t register_embree_ghost_dynamic_op::operator()(const TriangleMesh &mesh) const
{
    RTCGeometry rtc_geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
    assert(time_steps == 2); // we only support old and new frame for now
    rtcSetGeometryTimeStepCount(rtc_geom, 2);

    /// for the old frame flip the two meshes
    {
        const TriangleMesh &mesh_old = (const TriangleMesh &)ghost_dynamic.old_shape;
        Vector4f *positions = (Vector4f *)rtcSetNewGeometryBuffer(
            rtc_geom, RTC_BUFFER_TYPE_VERTEX, 1, RTC_FORMAT_FLOAT3,
            sizeof(Vector4f), mesh_old.positions.size());
        Vector3i *triangles = (Vector3i *)rtcSetNewGeometryBuffer(
            rtc_geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
            sizeof(Vector3i), mesh_old.indices.size());
        for (int i = 0; i < (int)mesh_old.positions.size(); i++)
        {
            Vector3 position = mesh_old.positions[i];
            positions[i] = Vector4f{(float)position[0], (float)position[1], (float)position[2], 0.f};
        }
        for (int i = 0; i < (int)mesh_old.indices.size(); i++)
        {
            triangles[i] = mesh_old.indices[i];
        }

    }
    /// for the current frame
    { // share the indices

        Vector4f *positions = (Vector4f *)rtcSetNewGeometryBuffer(
            rtc_geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
            sizeof(Vector4f), mesh.positions.size()); //TODO attention here is 0

        for (int i = 0; i < (int)mesh.positions.size(); i++)
        {
            Vector3 position = mesh.positions[i];
            positions[i] = Vector4f{(float)position[0], (float)position[1], (float)position[2], 0.f};
        }
    }

    ///////////////////////////
    // set intersection filter for the mesh
    {
        rtcSetGeometryIntersectFilterFunction(rtc_geom, intersection_filter);
        rtcSetGeometryOccludedFilterFunction(rtc_geom, occlusion_filter);
    }
    uint32_t geomID = rtcAttachGeometry(scene, rtc_geom);
    rtcAttachGeometryByID(ghost_scene, rtc_geom, geomID);
    rtcCommitGeometry(rtc_geom);
    rtcReleaseGeometry(rtc_geom);
    return geomID;
}

/////////////////////////// ghost_dynamic sphere//////////////
void ghost_dynamic_sphere_bounds_func(const struct RTCBoundsFunctionArguments *args)
{
}

void ghost_dynamic_sphere_intersect_func(const RTCIntersectFunctionNArguments *args)
{
    // TODO
}

void ghost_dynamic_sphere_occluded_func(const RTCOccludedFunctionNArguments *args)
{
    // TODO
}
uint32_t register_embree_ghost_dynamic_op::operator()(const Sphere &sphere) const
{
    RTCGeometry rtc_geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
    uint32_t geomID = rtcAttachGeometry(scene, rtc_geom);
    rtcAttachGeometryByID(ghost_scene, rtc_geom, geomID);

    rtcSetGeometryUserPrimitiveCount(rtc_geom, 1);
    rtcSetGeometryTimeStepCount(rtc_geom, 2);
    rtcSetGeometryUserData(rtc_geom, (void *)&ghost_dynamic);
    rtcSetGeometryBoundsFunction(rtc_geom, ghost_dynamic_sphere_bounds_func, nullptr);
    rtcSetGeometryIntersectFunction(rtc_geom, ghost_dynamic_sphere_intersect_func);
    rtcSetGeometryOccludedFunction(rtc_geom, ghost_dynamic_sphere_occluded_func);
    rtcCommitGeometry(rtc_geom);
    rtcReleaseGeometry(rtc_geom);
    return geomID;
}
