def format_components_info(components: list) -> str:
    if not components:
        return ""

    sorted_components = sorted(components, key=lambda x: x.get('order', 0))
    
    lines = []
    for comp in sorted_components:
        label = comp.get('label', '')
        order = comp.get('order', '')
        description = comp.get('description', '')
        lines.append(f"   组件 {label} (顺序: {order}): {description}")
    
    return '\n'.join(lines)


def format_example_json(components: list) -> str:
    if not components:
        return "{}"

    sorted_components = sorted(components, key=lambda x: x.get('order', 0))
    
    lines = ["{"]
    for i, comp in enumerate(sorted_components):
        label = comp.get('label', '')
        comma = "," if i < len(sorted_components) - 1 else ""
        lines.append(f'    "{label}": "具体内容"{comma}')
    lines.append("}")
    
    return '\n'.join(lines)