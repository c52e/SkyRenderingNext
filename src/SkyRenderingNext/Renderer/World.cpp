#include "World.h"

#include "Core/Common.h"
#include "Core/Gui.h"
#include "RenderManager.h"

static Entity gEditingEntity = entt::null;

World &GetWorld()
{
    static World s_World;
    return s_World;
}

entt::resource_cache<GltfModel, GltfModelLoader> &GetModelCache()
{
    static entt::resource_cache<GltfModel, GltfModelLoader> s_ModelCache{};
    return s_ModelCache;
}

void AddModel(GraphicsCommandList *cmd, const std::filesystem::path &path, const Math::float4x4 &matrix)
{
    auto id = entt::hashed_string(path.string().c_str());
    auto model = GetModelCache().load(id, cmd, path).first->second;

    auto &world = GetWorld();
    auto entity = world.create();
    world.emplace<Transform>(entity, matrix);
    world.emplace<ModelResourceHandle>(entity, model.handle());
    world.emplace<EntityName>(entity, path.stem().string());
    world.emplace<ResourcePath>(entity, path);
    gEditingEntity = entity;
}

static ImGuizmo::OPERATION mCurrentGizmoOperation(ImGuizmo::TRANSLATE);
static ImGuizmo::MODE mCurrentGizmoMode(ImGuizmo::WORLD);

// https://github.com/CedricGuillemet/ImGuizmo
bool EditTransform(glm::mat4 &matrix)
{
    bool changed = false;
    if (ImGui::RadioButton("Translate", mCurrentGizmoOperation == ImGuizmo::TRANSLATE))
        mCurrentGizmoOperation = ImGuizmo::TRANSLATE;
    ImGui::SameLine();
    if (ImGui::RadioButton("Rotate", mCurrentGizmoOperation == ImGuizmo::ROTATE))
        mCurrentGizmoOperation = ImGuizmo::ROTATE;
    ImGui::SameLine();
    if (ImGui::RadioButton("Scale", mCurrentGizmoOperation == ImGuizmo::SCALE))
        mCurrentGizmoOperation = ImGuizmo::SCALE;
    float matrixTranslation[3], matrixRotation[3], matrixScale[3];
    ImGuizmo::DecomposeMatrixToComponents(glm::value_ptr(matrix), matrixTranslation, matrixRotation, matrixScale);
    changed |= ImGui::InputFloat3("Tr", matrixTranslation);
    changed |= ImGui::InputFloat3("Rt", matrixRotation);
    changed |= ImGui::InputFloat3("Sc", matrixScale);
    ImGuizmo::RecomposeMatrixFromComponents(matrixTranslation, matrixRotation, matrixScale, glm::value_ptr(matrix));

    if (mCurrentGizmoOperation != ImGuizmo::SCALE)
    {
        if (ImGui::RadioButton("Local", mCurrentGizmoMode == ImGuizmo::LOCAL))
            mCurrentGizmoMode = ImGuizmo::LOCAL;
        ImGui::SameLine();
        if (ImGui::RadioButton("World", mCurrentGizmoMode == ImGuizmo::WORLD))
            mCurrentGizmoMode = ImGuizmo::WORLD;
    }

    return changed;
}

bool EditWorld()
{
    auto &world = GetWorld();
    if (!world.valid(gEditingEntity))
    {
        gEditingEntity = entt::null;
    }

    bool sceneDirty = false;
    ImGui::Begin("Entities");
    for (const auto &[entity, _] : world.view<Transform>().each())
    {
        // https://github.com/ocornut/imgui/issues/581
        ImGuiTreeNodeFlags node_flags = 0;
        node_flags |= (gEditingEntity == entity ? ImGuiTreeNodeFlags_Selected : 0);
        node_flags |= false ? ImGuiTreeNodeFlags_OpenOnArrow : ImGuiTreeNodeFlags_Leaf;
        auto name = world.try_get<EntityName>(entity);
        bool opened = ImGui::TreeNodeEx(reinterpret_cast<void *>(static_cast<intptr_t>(entity)), node_flags,
                                        name ? name->value.c_str() : "");
        if (ImGui::IsItemClicked())
        {
            gEditingEntity = entity;
        }
        else if (ImGui::IsItemClicked(ImGuiMouseButton_Right) && gEditingEntity == entity)
        {
            gEditingEntity = entt::null;
        }
        else if (ImGui::IsItemClicked(ImGuiMouseButton_Middle) && gEditingEntity == entity)
        {
            sceneDirty = true;
            world.destroy(entity);
            gEditingEntity = entt::null;
        }
        if (opened)
        {
            ImGui::TreePop();
        }
    }
    ImGui::End();

    ImGui::Begin("Detail");
    ImGuizmo::BeginFrame();
    if (gEditingEntity != entt::null)
    {
        auto &transform = world.get<Transform>(gEditingEntity);
        sceneDirty |= EditTransform(transform.m);
        auto windowRec = GetRenderManager().GetWindowAbsoluteRect();
        ImGuizmo::SetRect(windowRec.x, windowRec.y, windowRec.z - windowRec.x, windowRec.w - windowRec.y);
        auto &cam = GetRenderManager().MainCamera();
        sceneDirty |= ImGuizmo::Manipulate(
            Math::value_ptr(cam.ViewMatrix()), Math::value_ptr(cam.ProjectionFarClip(cam.m_zNear * 1e4f)),
            mCurrentGizmoOperation, mCurrentGizmoMode, Math::value_ptr(transform.m), NULL, NULL);
    }
    ImGui::End();
    return sceneDirty;
}
