import yaml
from pathlib import Path

# Define constants here
SCRIPT_DIR = Path(__file__).resolve().parent

WORKFLOW_DIR = SCRIPT_DIR / 'workflow'
DEFAULT_CONFIG = WORKFLOW_DIR / 'config' / 'config.yaml'
USER_CONFIG_DIR = Path.home() / ".config" / "edshat"
USER_CONFIG_FILE = USER_CONFIG_DIR / "config.yaml"


def _convert_value(value_str):
    val_lower = value_str.lower()
    if val_lower == 'true':
        return True
    if val_lower == 'false':
        return False
    try:
        return int(value_str)
    except ValueError:
        pass
    try:
        return float(value_str)
    except ValueError:
        pass
    if val_lower in ['none', 'null']:
        return None
    return value_str

def _convert_type(value_str):
    """Helper function to intelligently convert a string value to bool, int, or float."""
    if not isinstance(value_str, str):
        return value_str
        
    val_lower = value_str.lower()

    if val_lower.startswith('['):
        value_str = value_str.strip('[]')
        vals = []
        for val in value_str.split(','):
            vals.append(_convert_value(val.strip()))
        
        return vals
    return _convert_value(value_str)


def change_config(config, key, value):
    key_parts = key.split(';')
    d = config
    for part in key_parts[:-1]:
        # setdefault is perfect here: it gets the key or creates a new dict if it's missing
        d = d.setdefault(part, {})
    
    # Set the final key to the converted value
    d[key_parts[-1]] = _convert_type(value)
    

def set_config(args):
    """Sets a configuration variable and saves it to the user config file."""
    # Ensure the user config directory exists
    USER_CONFIG_DIR.mkdir(parents=True, exist_ok=True)

    # Load the existing user config, or start with an empty dict if it doesn't exist
    if USER_CONFIG_FILE.exists():
        with open(USER_CONFIG_FILE, 'r') as f:
            config = yaml.safe_load(f) or {}
    else:
        config = {}

    # Traverse the dictionary using the semicolon-separated key
    key_parts = args.variable.split(';')
    d = config
    for part in key_parts[:-1]:
        # setdefault is perfect here: it gets the key or creates a new dict if it's missing
        d = d.setdefault(part, {})
    
    # Set the final key to the converted value
    d[key_parts[-1]] = _convert_type(args.value)

    # Write the updated dictionary back to the user config file
    with open(USER_CONFIG_FILE, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    
    print(f"Configuration updated! '{args.variable}' set in {USER_CONFIG_FILE}")

def deep_update(d, u):
    """Recursively update a dictionary."""
    for k, v in u.items():
        if isinstance(v, dict):
            d[k] = deep_update(d.get(k, {}), v)
        else:
            d[k] = v
    return d

def get_config(args):
    """Loads configurations from default and user files."""
    # 1. Load the default configuration from within the package
    with open(DEFAULT_CONFIG, 'r') as f:
        config = yaml.safe_load(f)

    # 2. Load user's custom config if it exists and merge it
    if USER_CONFIG_FILE.exists():
        print(f"Loading user configuration from: {USER_CONFIG_FILE}")
        with open(USER_CONFIG_FILE, 'r') as f:
            user_config = yaml.safe_load(f)
        config = deep_update(config, user_config)
    
    # 3. Load configfile / config values from arguements
    if args.config or args.config_file:
        if args.config_file:
            with open(args.config_file, 'r') as f:
                args_config = yaml.safe_load(f)
            config = deep_update(config, args_config)
        if args.config:
            for item in args.config:
                key, value = item.split('=')
                change_config(config, key, value)

    return config
