#pragma once

#include "CommandList.h"
#include "Core/Common.h"

template <class ResourceType> class ResourceBarrierBase
{
  public:
    void BarrierTransition(GraphicsCommandList *cmd, D3D12_RESOURCE_STATES stateBefore,
                           D3D12_RESOURCE_STATES stateAfter, bool flush = true)
    {
        auto barrier = CD3DX12_RESOURCE_BARRIER::Transition(static_cast<ResourceType *>(this)->GetResource(),
                                                            stateBefore, stateAfter);
        cmd->ResourceBarrier(barrier, flush);
        m_ResourceState = stateAfter;
    }

    void BarrierTransitionTo(GraphicsCommandList *cmd, D3D12_RESOURCE_STATES stateAfter, bool flush = true)
    {
        BarrierTransition(cmd, m_ResourceState, stateAfter, flush);
    }

    void BarrierTransitionToCopySrc(GraphicsCommandList *cmd, bool flush = true)
    {
        BarrierTransitionTo(cmd, D3D12_RESOURCE_STATE_COPY_SOURCE, flush);
    }

    void BarrierTransitionToSrv(GraphicsCommandList *cmd, bool flush = true)
    {
        BarrierTransitionTo(
            cmd, D3D12_RESOURCE_STATE_NON_PIXEL_SHADER_RESOURCE | D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE, flush);
    }

    void BarrierTransitionToUav(GraphicsCommandList *cmd, bool flush = true)
    {
        BarrierTransitionTo(cmd, D3D12_RESOURCE_STATE_UNORDERED_ACCESS, flush);
    }

    void BarrierUAV(GraphicsCommandList *cmd, bool flush = true)
    {
        auto barrier = CD3DX12_RESOURCE_BARRIER::UAV(static_cast<ResourceType *>(this)->GetResource());
        cmd->ResourceBarrier(barrier, flush);
    }

  protected:
    ResourceBarrierBase(D3D12_RESOURCE_STATES state = D3D12_RESOURCE_STATE_COMMON) : m_ResourceState(state)
    {
    }

    D3D12_RESOURCE_STATES m_ResourceState;
};
