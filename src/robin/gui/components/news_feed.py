"""
News feed manager for ROBIN application.

This module handles fetching and caching news from the ROBIN news API.
"""

import json
import logging
from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import Any, Dict, List, Optional

import requests
from nicegui import run, ui

from robin.gui.theme import client_timer, stop_timer, ui_element_exists


@dataclass
class NewsRecord:
    id: str
    headline: str
    start_date: str
    end_date: str
    content: str
    link: Optional[str] = None
    image_url: Optional[str] = None
    message_type: str = "news"

    @classmethod
    def from_dict(cls, item: Dict[str, Any]) -> "NewsRecord":
        """Create a NewsRecord instance from a dictionary."""
        raw_image = item.get("image_url")
        if isinstance(raw_image, str):
            raw_image = raw_image.strip() or None
        return cls(
            id=item["id"],
            headline=item["headline"],
            start_date=item["start_date"],
            end_date=item["end_date"],
            content=item["content"],
            link=item.get("link"),
            image_url=raw_image,
            message_type=item.get("message_type") or "news",
        )


class NewsFeed:
    """Manages fetching and displaying news feed content."""

    def __init__(self):
        self.api_url = (
            "https://t373wby715.execute-api.eu-west-1.amazonaws.com/prod/news"
        )
        self.news_items: List[NewsRecord] = []
        self.last_update = None
        self.update_interval = timedelta(hours=6)
        self.is_available = True
        self._news_container = None
        self.headers = {
            "Content-Type": "application/json",
            "User-Agent": "ROBIN-Client/1.0",
            "Accept": "application/json",
        }
        self._update_timer = None

    async def fetch_news(self) -> bool:
        """
        Fetch news from the API.

        Returns:
            bool: True if successful, False otherwise
        """
        try:
            logging.debug("Fetching news feed...")
            response = await run.io_bound(
                requests.get,
                self.api_url,
                headers=self.headers,
                timeout=5,  # 5 second timeout
            )

            if response.status_code == 200:
                # Parse the response which has a nested structure
                data = response.json()
                # The body is a JSON string that needs to be parsed
                if isinstance(data.get("body"), str):
                    body_data = json.loads(data["body"])
                else:
                    body_data = data.get("body", {})

                # News items are in the 'news' array
                news_items = body_data.get("news", [])
                if news_items:
                    self.news_items = [
                        NewsRecord.from_dict(item) for item in news_items
                    ]
                    self.last_update = datetime.now()
                    self.is_available = True
                    logging.info(
                        f"Successfully fetched {len(self.news_items)} news items"
                    )
                    return True
                else:
                    logging.warning("No news items found in response")
                    self.is_available = True  # API is available, just no news
                    self.news_items = []
                    return True
            elif response.status_code == 403:
                logging.warning(
                    "Authentication failed when fetching news feed - please check API access"
                )
                self.is_available = False
                return False
            elif response.status_code == 404:
                logging.warning("News feed endpoint not found")
                self.is_available = False
                return False
            else:
                logging.warning(
                    f"Failed to fetch news. Status: {response.status_code}, Response: {response.text}"
                )
                self.is_available = False
                return False

        except requests.exceptions.Timeout:
            logging.warning("Timeout while fetching news feed")
            self.is_available = False
            return False
        except requests.exceptions.ConnectionError:
            logging.warning(
                "Network error while fetching news feed - check your internet connection"
            )
            self.is_available = False
            return False
        except json.JSONDecodeError as e:
            logging.error(f"Error parsing news feed response: {str(e)}")
            self.is_available = False
            return False
        except Exception as e:
            logging.error(f"Error fetching news: {str(e)}")
            self.is_available = False
            return False

    def should_update(self) -> bool:
        """Check if it's time to update the news feed."""
        if not self.last_update:
            return True
        return datetime.now() - self.last_update > self.update_interval

    def create_news_element(self) -> None:
        """Create the news feed UI element."""
        try:
            # Header section with title and last update
            with ui.row().classes("w-full"):
                with ui.row().classes("items-center gap-2"):
                    ui.icon("feed", color="primary").classes("text-xl")
                    ui.label("ROBIN News").classes(
                        "text-headline-medium text-slate-900 dark:text-slate-50"
                    )
                if self.last_update:
                    with ui.row().classes("items-center gap-1"):
                        ui.icon("update", color="gray").classes("text-sm")
                        ui.label(
                            f'Updated {self.last_update.strftime("%H:%M %d/%m/%Y")}'
                        ).classes("text-body-small text-slate-600 dark:text-slate-400")

            # Scrollable news container with elegant styling
            with ui.scroll_area().classes("w-full h-96"):
                self._news_container = ui.column().classes("w-full gap-3 p-3")
                self._update_news_display()

        except Exception as e:
            logging.error(f"Error creating news element: {str(e)}")
            with ui.row().classes("items-center gap-2 p-4"):
                ui.icon("error", color="negative")
                ui.label("Error Loading News").classes(
                    "text-xl font-bold text-negative"
                )

    def _news_container_alive(self) -> bool:
        """Return True if the news container still belongs to a live client."""
        return bool(self._news_container) and ui_element_exists(self._news_container)

    def _update_news_display(self) -> None:
        """Update the news display with current items."""
        # Timer callbacks can outlive the page (navigation / disconnect). Never
        # touch UI after the container's client has been deleted.
        if not self._news_container_alive():
            return

        try:
            self._news_container.clear()

            with self._news_container:
                if not self.is_available:
                    with ui.card().classes(
                        "w-full rounded-lg p-4 border border-[color:var(--md-error)]/30 bg-[color:var(--md-error-container)]"
                    ):
                        with ui.row().classes("items-center gap-2"):
                            ui.icon("wifi_off", color="negative")
                            ui.label(
                                "News feed unavailable - please check your network connection"
                            ).classes(
                                "text-body-medium text-[color:var(--md-on-error-container)]"
                            )
                    return

                if not self.news_items:
                    with ui.card().classes(
                        "w-full rounded-lg p-4 bg-slate-50 dark:bg-slate-900/40"
                    ):
                        with ui.row().classes("items-center gap-2"):
                            ui.icon("inbox", color="gray")
                            ui.label("No news items available").classes(
                                "text-body-medium text-slate-600 dark:text-slate-400"
                            )
                    return

                for item in self.news_items:
                    with ui.card().classes(
                        "w-full transition-shadow duration-200 hover:shadow-md"
                    ):
                        with ui.column().classes("w-full p-4 gap-3"):
                            # Headline with icon based on message type
                            with ui.row().classes("items-center gap-2"):
                                # Choose icon based on message type
                                icon_name = {
                                    "news": "newspaper",
                                    "update": "system_update",
                                    "alert": "warning",
                                }.get(item.message_type, "info")
                                ui.icon(icon_name, color="primary").classes("text-lg")
                                ui.label(item.headline).classes(
                                    "text-headline-small text-slate-900 dark:text-slate-50"
                                )

                            # Main content with proper spacing and formatting
                            with ui.column().classes(
                                "gap-3 pl-8"
                            ):  # Indent content under headline
                                ui.label(item.content).classes(
                                    "text-body-large text-slate-600 dark:text-slate-400 leading-relaxed"
                                )

                                # Image handling with better presentation
                                if item.image_url:
                                    try:
                                        # Verify image URL is valid
                                        if not item.image_url.startswith(
                                            ("http://", "https://")
                                        ):
                                            raise ValueError("Invalid image URL")

                                        # Create a full-width container with fixed dimensions
                                        with ui.card().classes(
                                            "w-full overflow-hidden rounded-lg my-1"
                                        ):
                                            with ui.row().classes(
                                                "justify-center items-center w-full"
                                            ):
                                                # Quasar QImg (ui.image) often collapses to ~0 height here until a
                                                # ratio is set; native <img> gets a real src on first paint and
                                                # sizes from the loaded bitmap like a normal browser image.
                                                ui.element("img").classes(
                                                    "block w-full max-w-md h-auto rounded-lg object-contain"
                                                ).props(
                                                    f"src={json.dumps(item.image_url)} "
                                                    f"alt={json.dumps(item.headline[:200])} "
                                                    "loading=lazy referrerpolicy=no-referrer"
                                                )
                                    except Exception as e:
                                        logging.warning(
                                            f"Error displaying image: {str(e)}"
                                        )
                                        with ui.row().classes(
                                            "items-center justify-center gap-2 p-2"
                                        ):
                                            ui.icon("broken_image", color="gray")
                                            ui.label("Image unavailable").classes(
                                                "text-gray-500 text-sm"
                                            )

                                # Link handling with improved styling
                                if item.link:
                                    with ui.row().classes("w-full justify-end mt-2"):
                                        with ui.link(
                                            target=item.link, new_tab=True
                                        ).classes(
                                            "flex items-center gap-1 text-primary hover:text-primary-dark transition-colors"
                                        ):
                                            ui.label("Read more")
                                            ui.icon("arrow_forward", size="text-sm")

                            # Add date information if available
                            if hasattr(item, "start_date"):
                                with ui.row().classes(
                                    "w-full justify-end mt-2 pt-2 border-t border-[color:var(--md-outline)]"
                                ):
                                    ui.label(f"Posted: {item.start_date}").classes(
                                        "text-label-small text-slate-500 dark:text-slate-400 font-mono"
                                    )

        except RuntimeError as e:
            # Client/slot torn down mid-update (navigation or disconnect).
            if "deleted" in str(e).lower() or "slot stack" in str(e).lower():
                return
            logging.error(f"Error updating news display: {e}")
            self._show_update_error()
        except Exception as e:
            logging.error(f"Error updating news display: {e}")
            self._show_update_error()

    def _show_update_error(self) -> None:
        """Show an in-feed error card only if the container is still usable."""
        if not self._news_container_alive():
            return
        try:
            with self._news_container:
                with ui.card().classes(
                    "w-full rounded-lg p-4 border border-[color:var(--md-error)]/30 bg-[color:var(--md-error-container)]"
                ):
                    with ui.row().classes("items-center gap-2"):
                        ui.icon("error", color="negative")
                        ui.label("Error updating news feed").classes(
                            "text-body-medium text-[color:var(--md-on-error-container)]"
                        )
        except RuntimeError as e:
            if "deleted" in str(e).lower() or "slot stack" in str(e).lower():
                return
            logging.debug(f"Could not show news update error UI: {e}")

    async def _check_and_update_news(self) -> None:
        """Check if update is needed and fetch news if necessary."""
        if self.should_update():
            if await self.fetch_news():
                self._update_news_display()

    def start_update_timer(self) -> None:
        """Start a client-scoped hourly update timer.

        Uses ``client_timer`` (app.timer + disconnect/delete cleanup) so navigation
        or tab close cannot leave a callback updating a deleted client. UI updates
        are also guarded by ``_news_container_alive``.
        """
        if self._update_timer is not None:
            return

        self._update_timer = client_timer(
            3600,
            self._check_and_update_news,
            immediate=True,
            active=True,
        )
        logging.info("News feed update timer initialized")

    def stop_update_timer(self) -> None:
        """Stop the periodic update timer if running."""
        stop_timer(self._update_timer)
        self._update_timer = None
