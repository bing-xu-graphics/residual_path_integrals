
struct register_embree_dynamic_op
{
    uint32_t operator()(const Sphere &sphere) const;
    uint32_t operator()(const TriangleMesh &mesh) const;
    const RTCDevice &device;
    const RTCScene &scene;
    const int time_steps;
    Dynamic &dynamic;
};

uint32_t register_embree_dynamic_op::operator()(const TriangleMesh &mesh) const
{
    RTCGeometry rtc_geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
    assert(time_steps == 2); // we only support old and new frame for now
    rtcSetGeometryTimeStepCount(rtc_geom, 2);
    
    /// for the old frame
    {
        const TriangleMesh &mesh_old = (const TriangleMesh &)dynamic.old_shape;
        Vector4f *positions = (Vector4f *)rtcSetNewGeometryBuffer(
            rtc_geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
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
            rtc_geom, RTC_BUFFER_TYPE_VERTEX, 1, RTC_FORMAT_FLOAT3,
            sizeof(Vector4f), mesh.positions.size());

        for (int i = 0; i < (int)mesh.positions.size(); i++)
        {
            Vector3 position = mesh.positions[i];
            positions[i] = Vector4f{(float)position[0], (float)position[1], (float)position[2], 0.f};
        }
    }
    uint32_t geomID = rtcAttachGeometry(scene, rtc_geom);
    rtcCommitGeometry(rtc_geom);
    rtcReleaseGeometry(rtc_geom);
    return geomID;
}

/////////////////////////// dynamic sphere//////////////
void dynamic_sphere_bounds_func(const struct RTCBoundsFunctionArguments *args)
{
    Dynamic *dynamic = (Dynamic *)args->geometryUserPtr;
    RTCBounds *bounds_o = args->bounds_o;
    const unsigned int time = args->timeStep;
    if (time > 0)
    {
        const Sphere *sphere = (const Sphere *)&dynamic->base_shape;
        bounds_o->lower_x = sphere->position.x - sphere->radius;
        bounds_o->lower_y = sphere->position.y - sphere->radius;
        bounds_o->lower_z = sphere->position.z - sphere->radius;
        bounds_o->upper_x = sphere->position.x + sphere->radius;
        bounds_o->upper_y = sphere->position.y + sphere->radius;
        bounds_o->upper_z = sphere->position.z + sphere->radius;
    }
    else
    {
        const Sphere *sphere = (const Sphere *)&dynamic->old_shape;
        bounds_o->lower_x = sphere->position.x - sphere->radius;
        bounds_o->lower_y = sphere->position.y - sphere->radius;
        bounds_o->lower_z = sphere->position.z - sphere->radius;
        bounds_o->upper_x = sphere->position.x + sphere->radius;
        bounds_o->upper_y = sphere->position.y + sphere->radius;
        bounds_o->upper_z = sphere->position.z + sphere->radius;
    }
}

void dynamic_sphere_intersect_func(const RTCIntersectFunctionNArguments *args)
{
    // TODO
}

void dynamic_sphere_occluded_func(const RTCOccludedFunctionNArguments *args)
{
    // TODO
}
uint32_t register_embree_dynamic_op::operator()(const Sphere &sphere) const
{
    RTCGeometry rtc_geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
    uint32_t geomID = rtcAttachGeometry(scene, rtc_geom);
    rtcSetGeometryUserPrimitiveCount(rtc_geom, 1);
    rtcSetGeometryTimeStepCount(rtc_geom, 2);
    rtcSetGeometryUserData(rtc_geom, (void *)&dynamic);
    rtcSetGeometryBoundsFunction(rtc_geom, dynamic_sphere_bounds_func, nullptr);
    rtcSetGeometryIntersectFunction(rtc_geom, dynamic_sphere_intersect_func);
    rtcSetGeometryOccludedFunction(rtc_geom, dynamic_sphere_occluded_func);
    rtcCommitGeometry(rtc_geom);
    rtcReleaseGeometry(rtc_geom);
    return geomID;
}
