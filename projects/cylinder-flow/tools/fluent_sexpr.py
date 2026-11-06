"""Minimal data-only reader for Fluent's saved Scheme lists (never executes Scheme)."""
import json
import re

TOKENS = re.compile(r'\s+|;[^\n]*|\(|\)|"(?:\\.|[^"\\])*"|[^\s()]+')


def read_all(text):
    stack = [[]]
    for match in TOKENS.finditer(text):
        token = match.group()
        if token.isspace() or token.startswith(';'):
            continue
        if token == '(':
            value = []
            stack[-1].append(value)
            stack.append(value)
        elif token == ')':
            if len(stack) == 1:
                raise ValueError('unbalanced closing parenthesis')
            stack.pop()
        elif token.startswith('"'):
            try:
                stack[-1].append(json.loads(token))
            except json.JSONDecodeError:
                stack[-1].append(token[1:-1])
        elif token in ('#t', '#f'):
            stack[-1].append(token == '#t')
        else:
            try:
                value = int(token)
            except ValueError:
                try:
                    value = float(token)
                except ValueError:
                    value = token
            stack[-1].append(value)
    if len(stack) != 1:
        raise ValueError('unbalanced opening parenthesis')
    return stack[0]


def pairs(items):
    return {item[0]: item[1:] for item in items if isinstance(item, list) and item and isinstance(item[0], str)}
