#include "CommandList.h"

GraphicsCommandList::GraphicsCommandList(ID3D12Device *device, ID3D12CommandAllocator* alloc)
{
    device->CreateCommandList(0, D3D12_COMMAND_LIST_TYPE_DIRECT, alloc, nullptr, IID_PPV_ARGS(&m_CmdList));
}
